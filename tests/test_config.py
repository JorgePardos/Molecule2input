"""Profile loading, validation and CLI merging."""

from __future__ import annotations

import pytest

from m2i import config


def test_all_builtin_profiles_are_valid():
    profiles = config.available_profiles()
    assert profiles
    for name in profiles:
        config.load_profile(name).validate()


def test_unknown_profile_lists_what_exists():
    with pytest.raises(config.ProfileError) as excinfo:
        config.load_profile("does_not_exist")
    assert "gaussian_opt_freq" in str(excinfo.value)


def test_unknown_key_is_rejected():
    with pytest.raises(config.ProfileError):
        config.JobProfile.from_dict({"methodd": "b3lyp"})


def test_unknown_program_is_rejected():
    with pytest.raises(config.ProfileError):
        config.JobProfile.from_dict({"program": "gamess"})


def test_unknown_dispersion_keyword_is_rejected():
    with pytest.raises(config.ProfileError):
        config.JobProfile.from_dict({"program": "gaussian", "dispersion": "d3bj"})
    # the same keyword is correct for ORCA
    config.JobProfile.from_dict({"program": "orca", "dispersion": "d3bj"})


def test_unknown_solvent_model_is_rejected():
    with pytest.raises(config.ProfileError):
        config.JobProfile.from_dict(
            {"program": "gaussian", "solvent": {"model": "cosmo", "name": "water"}}
        )


def test_cli_flags_win_over_the_profile():
    profile = config.load_profile("gaussian_opt_freq")
    merged = profile.merged_with(basis="def2TZVP", method=None)
    assert merged.basis == "def2TZVP"
    assert merged.method == profile.method  # None means "not overridden"


def test_switching_program_translates_the_dispersion_keyword():
    """--program orca on a Gaussian profile must not die on 'gd3bj'."""
    gaussian = config.load_profile("gaussian_opt_freq")
    assert gaussian.dispersion == "gd3bj"
    assert gaussian.merged_with(program="orca").dispersion == "d3bj"

    orca = config.load_profile("orca_opt_freq")
    assert orca.merged_with(program="gaussian").dispersion == "gd3bj"


def test_dispersion_is_dropped_for_geometry_only_output():
    profile = config.load_profile("gaussian_opt_freq").merged_with(program="xyz")
    assert profile.dispersion is None


def test_dispersion_without_an_equivalent_is_reported_not_silently_dropped():
    profile = config.load_profile("orca_opt_freq").merged_with(dispersion="d4")
    with pytest.raises(config.ProfileError) as excinfo:
        profile.merged_with(program="gaussian")
    assert "d4" in str(excinfo.value)


def test_explicit_dispersion_override_is_not_translated():
    profile = config.load_profile("gaussian_opt_freq")
    assert profile.merged_with(program="orca", dispersion="d4").dispersion == "d4"


def test_local_profiles_shadow_builtin_ones(tmp_path, monkeypatch):
    directory = tmp_path / "profiles"
    directory.mkdir()
    (directory / "gaussian_opt_freq.yaml").write_text(
        "program: gaussian\nmethod: wb97xd\nbasis: def2SVP\njobs: [opt]\n"
    )
    monkeypatch.setenv("M2I_PROFILE_PATH", str(directory))
    assert config.load_profile("gaussian_opt_freq").method == "wb97xd"


@pytest.mark.parametrize(
    "text, megabytes",
    [("8GB", 8192), ("8 GB", 8192), ("512MB", 512), ("2", 2048), ("1TB", 1048576)],
)
def test_memory_parsing(text, megabytes):
    assert config._mem_to_mb(text) == megabytes


def test_orca_maxcore_is_derived_from_memory_and_cores():
    resources = config.Resources(mem="8GB", nproc=8)
    assert resources.orca_maxcore() == 768
    assert config.Resources(mem="8GB", nproc=8, maxcore_mb=2000).orca_maxcore() == 2000
