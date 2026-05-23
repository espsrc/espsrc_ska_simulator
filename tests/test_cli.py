"""CLI configuration behavior."""

from skasim import cli
import skasim.pipeline
from pydantic import ValidationError
import pytest


def test_catalog_sets_named_catalog(monkeypatch):
    """--catalog is the canonical CLI spelling for built-in catalogs."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--catalog", "GLEAM"])

    assert captured[0].catalog == "GLEAM"


def test_numeric_catalog_cli_value_has_migration_message():
    """Numeric catalog CLI values fail with the 0.2 migration message."""
    with pytest.raises(ValidationError, match="Numeric catalog IDs were removed"):
        cli.main(["--catalog", "1"])


def test_flux_density_sets_generated_source_flux_jy(monkeypatch):
    """--flux-density defines generated source intensities."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(["--flux-density", "1.0", "5.0", "10.0"])

    assert captured[0].source_flux_jy == [1.0, 5.0, 10.0]


def test_generated_source_stokes_flags(monkeypatch):
    """--stokes-q/u/v define generated source polarization values."""
    captured = []
    monkeypatch.setattr(skasim.pipeline, "run", lambda config: captured.append(config))

    cli.main(
        [
            "--flux-density",
            "2.5",
            "3.5",
            "--stokes-q",
            "0.1",
            "0.2",
            "--stokes-u",
            "0.0",
            "0.3",
            "--stokes-v",
            "0.0",
            "-0.1",
        ]
    )

    assert captured[0].source_flux_jy == [2.5, 3.5]
    assert captured[0].stokes_q_jy == [0.1, 0.2]
    assert captured[0].stokes_u_jy == [0.0, 0.3]
    assert captured[0].stokes_v_jy == [0.0, -0.1]


def test_stokes_i_has_migration_message(capsys):
    """--stokes-i fails with a migration message for --flux-density."""
    with pytest.raises(SystemExit):
        cli.main(["--stokes-i", "2.5", "3.5"])

    captured = capsys.readouterr()
    assert "--flux-density" in captured.err


def test_legacy_i_flag_has_migration_message(capsys):
    """--I fails with a migration message for the clearer source intensity flags."""
    with pytest.raises(SystemExit):
        cli.main(["--I", "1.0"])

    captured = capsys.readouterr()
    assert "--flux-density" in captured.err


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


def test_cli_rejects_source_flux_jy_with_catalog():
    """CLI users cannot combine generated-source intensities with a catalog."""
    with pytest.raises(ValidationError, match="generated source mode"):
        cli.main(["--catalog", "MIGHTEE", "--flux-density", "1.0"])


def test_cli_help_has_single_canonical_surface(capsys):
    """Public help exposes canonical flags and hides removed aliases."""
    with pytest.raises(SystemExit):
        cli.main(["--help"])

    captured = capsys.readouterr()
    help_text = captured.out
    assert "--frequency-mhz" in help_text
    assert "--observation-time" in help_text
    assert "--output-dir" in help_text
    assert "--clean-iterations" in help_text
    assert "--flux-density" in help_text
    assert "--stokes-q" in help_text
    assert "--stokes-u" in help_text
    assert "--stokes-v" in help_text
    assert "--catalog" in help_text

    for removed in (
        "--stokes-i",
        "--json",
        "--json-fg",
        "--freq ",
        "--seconds",
        "--prefix",
        "--niter",
        "--Q",
        "--U",
        "--V",
    ):
        assert removed not in help_text
