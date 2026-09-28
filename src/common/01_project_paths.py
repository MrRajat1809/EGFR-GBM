"""Shared filesystem paths for the standalone analysis scripts."""
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def project_path(relative_path):
    """Return an absolute path without depending on the launch directory."""
    return str(PROJECT_ROOT / relative_path)


def use_working_directory(relative_path):
    """Keep relative outputs in their analysis folder, away from source code."""
    directory = PROJECT_ROOT / relative_path
    directory.mkdir(parents=True, exist_ok=True)
    os.chdir(directory)
