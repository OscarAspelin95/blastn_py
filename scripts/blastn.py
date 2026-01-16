import json
import tempfile
from multiprocessing import cpu_count
from pathlib import Path

import pandas as pd
from sh import blastn, makeblastdb

from scripts.parser import BlastCols, BlastConfig, get_best_hits, read_blast_tsv
from scripts.schema import BlastResult
from scripts.utils import _ensure_dir, _ensure_file


def get_blastdb(subject: Path, wd: Path) -> Path:
    db_base = subject.name

    db_dir = wd / "blastdb"
    _ensure_dir(db_dir)

    db_path = db_dir / db_base
    makeblastdb("-dbtype", "nucl", "-title", db_base, "-out", db_path, "-in", subject)

    return db_path


def run_blastn(query: Path, blast_db: Path, cfg: BlastConfig, wd: Path) -> Path:
    result_dir = wd / "blastn"
    _ensure_dir(result_dir)

    if (blast_tsv := result_dir / "blast_raw.tsv").is_file():
        return blast_tsv

    blast_cols = BlastCols.blast_str()
    blastn(
        "-query",
        query,
        "-db",
        blast_db,
        "-num_threads",
        cpu_count(),
        "-perc_identity",
        int(100 * cfg.frac_identity),
        "-outfmt",
        f"6 {blast_cols}",
        "-out",
        blast_tsv,
    )

    _ensure_file(blast_tsv)
    return blast_tsv


def write_results(blast_result: BlastResult, outdir: Path) -> Path:
    with (result_json := outdir / "blast_result.json").open("w") as f:
        json.dump(blast_result.model_dump(mode="json"), f, indent=4)

    _ensure_file(result_json)

    return result_json


def validate_results(best_hits_df: pd.DataFrame, cfg: BlastConfig) -> BlastResult:
    return BlastResult.model_validate(
        {"hits": best_hits_df.to_dict(orient="records"), "config": cfg}
    )


def nucleotide_blast(
    query: str, subject: str, cfg: BlastConfig, outdir: Path
) -> tuple[BlastResult, Path]:
    """ """

    with tempfile.TemporaryDirectory() as wd:
        wd = Path(wd)

        _ensure_file(query)
        _ensure_file(subject)

        # Actual BLASTing
        blast_db = get_blastdb(subject, wd)
        blast_tsv = run_blastn(query, blast_db, cfg, wd)
        blast_df = read_blast_tsv(blast_tsv, cfg)

        # Parse and validate results
        best_hits_df = get_best_hits(blast_df, margin=cfg.margin)
        blast_result = validate_results(best_hits_df, cfg)
        result_json = write_results(blast_result, outdir)

        return blast_result, result_json
