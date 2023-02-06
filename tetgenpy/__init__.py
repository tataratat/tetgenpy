import sys

from tetgenpy import _version, gentet, plc
from tetgenpy import tetgenpy_core as core
from tetgenpy.gentet import tetgen_exe, tetrahedralize
from tetgenpy.plc import PLC
from tetgenpy.tetgenpy_core import TetgenIO

__version__ = _version.__version__

__all__ = [
    "_version",
    "__version__",
    "plc",
    "core",
    "gentet",
    "tetrahedralize",
    "PLC",
    "tetgen_exe",
    "tetgen",
    "TetgenIO",
]


def _tetgen():
    """
    tetgen executable equivalent function.

    Parameters
    ----------
    None

    Returns
    -------
    exit_code: int
    """
    return tetgen_exe(sys.argv)
