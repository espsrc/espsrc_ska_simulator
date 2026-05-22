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


def test_flux_density_sets_generated_source_intensities(monkeypatch):
    """--flux-density defines generated source intensities."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--flux-density", "1.0", "5.0", "10.0"])

    assert captured[0].source_intensities == [1.0, 5.0, 10.0]


def test_stokes_i_alias_sets_generated_source_intensities(monkeypatch):
    """--stokes-i is an alias for generated source intensities."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--stokes-i", "2.5", "3.5"])

    assert captured[0].source_intensities == [2.5, 3.5]


def test_legacy_i_flag_has_migration_message(capsys):
    """--I fails with a migration message for the clearer source intensity flags."""
    with pytest.raises(SystemExit):
        cli.main(["--I", "1.0"])

    captured = capsys.readouterr()
    assert "--flux-density" in captured.err
    assert "--stokes-i" in captured.err
