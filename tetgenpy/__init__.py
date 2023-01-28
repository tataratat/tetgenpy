import sys

from tetgenpy import _version, plc
from tetgenpy import tetgenpy_core as core
from tetgenpy import tetrahedralize
from tetgenpy.plc import PLC

__version__ = _version.__version__

__all__ = [
    "_version",
    "__version__",
    "plc",
    "core",
    "tetrahedralize",
    "PLC",
    "tetgen_exe",
    "tetgen",
]


def tetgen_exe(argv):
    """
    Calls tetgen executable with argv as a list.
    This call performs the same routine as tetgen's main().
    The first element of argv is usually the executable name and tetgen also
    expects that. However, this can be a dummy.

    Parameters
    ----------
    argv: list

    Returns
    -------
    exit_code: int
      Will retun 0 if everything went well. -1 with RuntimeError without
      any informative error message.
    """
    try:
        # only list
        if not isinstance(argv, list):
            raise TypeError("argv should be a list.")

        if len(argv) == 0:
            raise ValueError("argv needs at least one element")

        return core.main(argv)

    except RuntimeError as err:
        # tetgen's termination/exit is caught as an error, which causes
        # messsy prints.
        # Informative messages would have been printed from tetgen by now.
        # Mute RuntimeError with default message.
        if str(err) == "Caught an unknown exception!":
            return -1
        else:
            raise err


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
