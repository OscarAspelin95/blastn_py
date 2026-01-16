import pytest

from scripts.schema import BlastConfig


def test_blast_config_default_margin():
    """Test that BlastConfig has a default margin of 0."""
    config = BlastConfig()
    assert config.margin == 0


def test_blast_config_valid_margin():
    """Test that BlastConfig accepts valid margin values."""
    config = BlastConfig(margin=5)
    assert config.margin == 5

    config = BlastConfig(margin=0)
    assert config.margin == 0

    config = BlastConfig(margin=100)
    assert config.margin == 100


def test_blast_config_invalid_margin():
    """Test that BlastConfig rejects negative margin values."""
    with pytest.raises(ValueError, match="is not >= 0"):
        BlastConfig(margin=-1)

    with pytest.raises(ValueError, match="is not >= 0"):
        BlastConfig(margin=-100)
