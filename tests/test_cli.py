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


def test_default_imager_is_oskar_dirty(monkeypatch):
    """The CLI defaults to the OSKAR dirty imager."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main([])

    assert captured[0].imaging.imager == "oskar-dirty"


def test_wsclean_imager_is_selected_explicitly(monkeypatch):
    """--imager wsclean selects WSClean image production."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--imager", "wsclean"])

    assert captured[0].imaging.imager == "wsclean"
    assert captured[0].imaging.algorithm == "wsclean_clean"


def test_cleaning_flag_has_migration_message(capsys):
    """--cleaning fails with a migration message for --imager wsclean."""
    with pytest.raises(SystemExit):
        cli.main(["--cleaning"])

    captured = capsys.readouterr()
    assert "--imager wsclean" in captured.err


def test_wsclean_command_cli_override(monkeypatch):
    """--wsclean-command stores the configured WSClean command."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))
    command = "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"

    cli.main(["--imager", "wsclean", "--wsclean-command", command])

    assert captured[0].imaging.wsclean_command == command
