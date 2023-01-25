from tetgenpy import core


def tetrahedralize(switches, tetgenio_in):
    out = core.TetgenIO()
    core.tetrahedralize(switches, tetgenio_in, out)
    return out
