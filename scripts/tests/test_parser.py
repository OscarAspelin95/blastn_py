import pytest

from scripts.parser import is_in_bounds


@pytest.mark.parametrize(
    "first,second,expected_in_bounds",
    [
        ((1, 1), (2, 2), False),
        ((1, 1), (1, 1), True),
        ((1, 5), (2, 5), True),
        ((1, 5), (1, 6), True),
        ((10, 20), (15, 25), True),
    ],
    ids=[
        "01. Non-overlapping",
        "02. Identical",
        "03. Subset",
        "04. Subset (should never happen).",
        "05. Overlapping",
    ],
)
def test_is_in_bounds(first: tuple[int, int], second: tuple[int, int], expected_in_bounds: bool):
    """Test is_in_bounds with default margin=0 (original behavior)."""
    start_1, end_1 = first
    start_2, end_2 = second
    assert is_in_bounds(start_1, end_1, start_2, end_2) == expected_in_bounds


@pytest.mark.parametrize(
    "first,second,margin,expected_in_bounds",
    [
        # margin=0 cases (same as original behavior)
        ((1, 1), (2, 2), 0, False),
        ((1, 1), (1, 1), 0, True),
        ((1, 10), (11, 20), 0, False),  # Adjacent, no overlap
        ((1, 10), (10, 20), 0, True),  # Touch at position 10 (overlap = 1)
        ((1, 10), (9, 15), 0, True),  # Overlap by 2 positions (9, 10)

        # margin=1 cases
        ((1, 10), (10, 20), 1, False),  # Overlap = 1, margin = 1, so separate
        ((1, 10), (9, 15), 1, True),  # Overlap = 2 > margin = 1, so overlapping
        ((1, 10), (8, 15), 1, True),  # Overlap = 3 > margin = 1, so overlapping
        ((1, 5), (5, 10), 1, False),  # Overlap = 1 = margin, so separate
        ((1, 1), (1, 1), 1, False),  # Overlap = 1 = margin, so separate

        # margin=2 cases (user's example)
        ((1, 10), (9, 15), 2, False),  # Overlap = 2 = margin, so separate
        ((1, 10), (8, 15), 2, True),  # Overlap = 3 > margin = 2, so overlapping
        ((1, 10), (10, 20), 2, False),  # Overlap = 1 < margin = 2, so separate
        ((1, 10), (11, 20), 2, False),  # No overlap

        # margin=3 cases
        ((1, 10), (8, 15), 3, False),  # Overlap = 3 = margin, so separate
        ((1, 10), (7, 15), 3, True),  # Overlap = 4 > margin = 3, so overlapping

        # Large margin cases
        ((1, 10), (2, 9), 10, False),  # Complete subset, overlap = 8 < margin = 10
        ((1, 100), (50, 150), 10, True),  # Large overlap = 51 > margin = 10
        ((1, 10), (5, 15), 10, False),  # Overlap = 6 < margin = 10

        # Edge cases with margin
        ((1, 1), (2, 2), 1, False),  # No overlap, any margin
        ((1, 5), (1, 5), 5, False),  # Identical, overlap = 5 = margin
        ((1, 5), (1, 5), 4, True),  # Identical, overlap = 5 > margin = 4
    ],
    ids=[
        "01. margin=0: adjacent no overlap",
        "02. margin=0: identical single position",
        "03. margin=0: adjacent intervals",
        "04. margin=0: touch at one position",
        "05. margin=0: overlap by 2",
        "06. margin=1: overlap=1 exactly",
        "07. margin=1: overlap=2 > margin",
        "08. margin=1: overlap=3 > margin",
        "09. margin=1: touch at edge",
        "10. margin=1: identical single position",
        "11. margin=2: user example (1,10) and (9,15)",
        "12. margin=2: overlap=3 > margin",
        "13. margin=2: overlap=1 < margin",
        "14. margin=2: adjacent no overlap",
        "15. margin=3: overlap=3 exactly",
        "16. margin=3: overlap=4 > margin",
        "17. margin=10: complete subset",
        "18. margin=10: large overlap",
        "19. margin=10: moderate overlap < margin",
        "20. edge: no overlap any margin",
        "21. edge: identical overlap = margin",
        "22. edge: identical overlap > margin",
    ],
)
def test_is_in_bounds_with_margin(
    first: tuple[int, int], second: tuple[int, int], margin: int, expected_in_bounds: bool
):
    """Test is_in_bounds with various margin values.

    The margin parameter allows hits that overlap by <= margin nucleotides
    to be considered separate locations.
    """
    start_1, end_1 = first
    start_2, end_2 = second
    assert is_in_bounds(start_1, end_1, start_2, end_2, margin) == expected_in_bounds
