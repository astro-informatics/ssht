"""SSHT module."""

from pathlib import Path

__doc__ = (Path(__file__).parent / "SSHT_Python_Documentation.md").read_text()

from pyssht.exceptions import ssht_spin_error, ssht_input_error  # type: ignore
from pyssht.cpyssht import *
