from pathlib import Path

import numpy as np
import pandas as pd

from scripts.schema import BlastCols, BlastConfig


def read_blast_tsv(blast_tsv: Path, cfg: BlastConfig) -> pd.DataFrame:
    """Read BLAST tsv and remove low quality hits."""
    df = pd.read_csv(blast_tsv, names=BlastCols.to_list(), sep="\t")

    sstart = np.minimum(df["sstart"], df["send"])
    send = np.maximum(df["sstart"], df["send"])

    return df.assign(
        frac_aligned=df["length"] / df[["slen", "qlen"]].min(axis=1),
        pident=df["pident"].div(100),
        sstart=sstart,
        send=send,
    ).query("frac_aligned >= @cfg.frac_aligned and pident >= @cfg.frac_identity")


def is_in_bounds(start_1: int, end_1: int, start_2: int, end_2: int, margin: int = 0) -> bool:
    """Check if two intervals overlap by more than the margin.

    Args:
        start_1: Start position of first interval
        end_1: End position of first interval
        start_2: Start position of second interval
        end_2: End position of second interval
        margin: Margin in nucleotides. Intervals overlapping by <= margin are considered separate.

    Returns:
        True if intervals overlap by more than margin, False otherwise.
    """
    # Check if intervals overlap at all
    if not (start_1 <= end_2 and start_2 <= end_1):
        return False

    # Calculate overlap size (number of positions that overlap)
    overlap = min(end_1, end_2) - max(start_1, start_2) + 1

    # Intervals are considered overlapping if overlap > margin
    return overlap > margin


def assign_hit_location(
    s: pd.Series, locations: list[tuple[int, int]], margin: int = 0
) -> int:
    for i, (location_start, location_end) in enumerate(locations, start=1):
        sstart, ssend = int(s["sstart"]), int(s["send"])

        if is_in_bounds(location_start, location_end, sstart, ssend, margin):
            return i

    raise ValueError("Unknown hit location")


def get_best_hit(s: pd.DataFrame) -> pd.Series:
    return s.sort_values(
        by=["gaps", "frac_aligned", "pident"], ascending=[True, False, False]
    ).iloc[0]


def get_best_hit_per_hit_location(s: pd.DataFrame, margin: int = 0) -> pd.DataFrame:
    if s.empty:
        return pd.DataFrame()

    s = s.sort_values(by=["sstart", "send"], ascending=[True, False])

    locations: list[tuple[int, int]] = []

    s_first = s.iloc[0]
    start, end = int(s_first["sstart"]), int(s_first["send"])

    for _, row in s.iterrows():
        current_start, current_end = int(row["sstart"]), int(row["send"])

        # We have a new hit location to check
        match is_in_bounds(start, end, current_start, current_end, margin):
            case True:
                start = min(start, current_start)
                end = max(end, current_end)
            case False:
                locations.append((start, end))
                start = current_start
                end = current_end

    locations.append((start, end))

    s["hit_location"] = s.apply(assign_hit_location, args=(locations, margin), axis=1)
    return s.groupby(by="hit_location").apply(get_best_hit, include_groups=False)


def get_best_hits(blast_df: pd.DataFrame, margin: int = 0) -> pd.DataFrame:
    return (
        blast_df.groupby(by=["sseqid", "sstrand", "sframe"])
        .apply(get_best_hit_per_hit_location, margin=margin, include_groups=True)
        .reset_index(drop=True)
    )
