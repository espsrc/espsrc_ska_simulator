"""CLI configuration behavior."""

from skasim import cli
import skasim.pipeline
from pydantic import ValidationError
import pytest


def test_catalog_alias_maps_to_named_catalogue(monkeypatch):
    """--catalog is a CLI spelling alias for --catalogue."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--catalog", "GLEAM"])

    assert captured[0].catalogue == "GLEAM"


def test_numeric_catalogue_cli_value_has_migration_message():
    """Numeric catalogue CLI values fail with the 0.2 migration message."""
    with pytest.raises(ValidationError, match="Numeric catalogue IDs were removed"):
        cli.main(["--catalogue", "1"])
