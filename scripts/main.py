import argparse
from pathlib import Path

from .blastn import nucleotide_blast
from .schema import BlastConfig
from .utils import _ensure_dir, _ensure_file
from yaspin import yaspin


def load_config(args) -> BlastConfig:
    cfg = BlastConfig()
    if args.min_identity:
        cfg.frac_identity = args.min_identity
    if args.min_aligned:
        cfg.frac_aligned = args.min_aligned
    if args.word_size:
        cfg.word_size = args.word_size
    return cfg


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-q", "--query", type=Path, required=True, help="")
    parser.add_argument("-s", "--subject", type=Path, required=True, help="")
    parser.add_argument("-o", "--outdir", type=Path, required=True)
    parser.add_argument("--min-identity", type=float, required=False)
    parser.add_argument("--min-aligned", type=float, required=False)
    parser.add_argument("--word-size", type=int, required=False)

    args = parser.parse_args()

    query = _ensure_file(args.query)
    subject = _ensure_file(args.subject)
    outdir = _ensure_dir(args.outdir)

    cfg = load_config(args)

    with yaspin(text="Running BLASTN...", color="green") as spinner:
        _, result_json = nucleotide_blast(query, subject, cfg, outdir)
        spinner.ok("✔")


if __name__ == "__main__":
    main()
