import argparse
import logging
from enum import Enum, unique
from pathlib import Path

import requests
from pydantic import BaseModel


@unique
class DataBaseType(Enum):
    AMR_FINDER_PLUS = "amr_finder_plus"

    @classmethod
    def as_list(cls) -> list[str]:
        return [db.value for db in cls]


class Database(BaseModel):
    url: str


class AmrFinderPlus(Database):
    url: str = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest/AMR_CDS.fa"

    def download(self, outdir: Path) -> Path:
        outdir.mkdir(exist_ok=True, parents=True)
        outfile = outdir / "amr_finder_plus.fasta"

        logging.info(f"Downloading {self.url}")

        response = requests.get(self.url)
        if response.status_code == 200:
            with outfile.open("wb") as f:
                f.write(response.content)
        else:
            raise Exception(f"Failed to download {self.url}")

        logging.info(f"Database file: {outfile}")
        return outfile


def get_database(database: DataBaseType) -> Database:
    match database:
        case DataBaseType.AMR_FINDER_PLUS:
            return AmrFinderPlus()
        case _:
            raise ValueError(f"Unknown database type: {database}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(description="Download databases")
    parser.add_argument(
        "-d",
        "--database",
        choices=DataBaseType.as_list(),
        help="Database to download",
        required=True,
    )
    parser.add_argument(
        "-o", "--outdir", type=Path, help="Output directory", default=Path.cwd()
    )
    args = parser.parse_args()

    database = get_database(DataBaseType(args.database))
    db_fasta = database.download(args.outdir)
