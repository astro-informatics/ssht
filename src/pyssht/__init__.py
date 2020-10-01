"""SSHT module."""

from pathlib import Path

__doc__ = (Path(__file__).parent / "SSHT_Python_Documentation.md").read_text()

from .pyssht import *
