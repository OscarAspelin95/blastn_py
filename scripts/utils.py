import re
from pathlib import Path

ACC_PAT = re.compile(
    r"^(?:[A-Z]{1,2}\d{5,6}|[A-Z]{2}_[A-Z]{2}\d{6,9}|[A-Z]{2}_\d{6,9}|[A-Z]{4}\d{8,9})(?:\.\d+)?$"
)


def valid_accession(accession: str) -> bool:
    return ACC_PAT.match(accession) is not None


def _ensure_file(f: str | Path) -> Path:
    if not (f_path := Path(f)).is_file():
        raise FileNotFoundError(f"File '{f}' not found.")
    return f_path


def _ensure_dir(d: str | Path) -> Path:
    if not (d_path := Path(d)).is_dir():
        d_path.mkdir(parents=True)
    return d_path
