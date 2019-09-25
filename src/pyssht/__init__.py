"""SSHT module."""

from pathlib import Path
__doc__ = (Path(__file__).parent / "SSHT_Python_documentation.md").read_text()

from .pyssht import *
