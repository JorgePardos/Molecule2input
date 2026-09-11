"""Make the package importable without installing it."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from m2i.types import IssueLog  # noqa: E402


@pytest.fixture
def log() -> IssueLog:
    return IssueLog()


@pytest.fixture(autouse=True)
def isolated_backend_home(monkeypatch, tmp_path_factory):
    """Never let the suite see the machine's real backend installations.

    Without this, whether a test passes depends on whether the developer
    happens to have run `m2i setup` -- and the failure appears in a test that
    says nothing about backends.
    """
    monkeypatch.setenv(
        "M2I_BACKEND_HOME", str(tmp_path_factory.mktemp("backend-home"))
    )


def codes(log: IssueLog) -> set[str]:
    """Issue codes present in a log, for concise assertions."""
    return {issue.code for issue in log}
