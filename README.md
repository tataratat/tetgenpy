# tetgenpy
[![main](https://github.com/tataratat/tetgenpy/actions/workflows/main.yml/badge.svg)](https://github.com/tataratat/tetgenpy/actions/workflows/main.yml)
[![PyPI version](https://badge.fury.io/py/tetgenpy.svg)](https://badge.fury.io/py/tetgenpy)

tetgenpy is a python wrapper for [Hang Si](https://www.wias-berlin.de/people/si/)'s [TetGen - A Quality Tetrahedral Mesh Generator and a 3D Delaunay Triangulator](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1).
It helps to prepare [various types of inputs - points, piecewise linear complexes (PLCs), and tetmesh -](https://wias-berlin.de/software/tetgen/1.5/doc/manual/manual027.png) for tetrahedron mesh generation based on simple python types, such as `list` and `numpy.ndarray`.

## Install
```bash
pip install tetgenpy
```
For current development version,
```bash
pip install git+https://github.com/tataratat/tetgenpy.git@main
```

## Quick Start
Following is an example for tetrahedralization of a unit cube defined as PLCs.
Alternatively, you could also use [tetgenpy.PLC](https://tataratat.github.io/tetgenpy/tetgenpy.html#module-tetgenpy.plc) class to prepare `TetgenIO`.
```python
import tetgenpy
import numpy as np

# tetrahedralize unit cube
# define points
points=[
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [0.0, 1.0, 1.0],
    [1.0, 1.0, 1.0],
]

# define facets - it can be list of polygons.
# here, they are hexa faces
facets = [
    [1, 0, 2, 3],
    [0, 1, 5, 4],
    [2, 0, 4, 6],
    [1, 3, 7, 5],
    [3, 2, 6, 7],
    [4, 5, 7, 6],
]

# prepare TetgenIO - input for tetgen
tetgen_in = tetgenpy.TetgenIO()

# set points, facets, and facet_markers.
# facet_markers can be useful for setting boundary conditions
tetgen_in.setup_plc(
    points=points,
    facets=facets,
    facet_markers=[[i] for i in range(1, len(facets) + 1)],
)

# tetgen's tetraheralize function with switches
tetgen_out = tetgenpy.tetrahedralize("qa.05", tetgen_in)

# unpack output
print(tetgen_out.points())
print(tetgen_out.tetrahedra())
print(tetgen_out.trifaces())
print(tetgen_out.trifacemarkers())
```
This package also provides access to tetgen's binary executable. Try:
```bash
$ tetgen -h
```

## Dependencies
#### c++
- [tetgen](https://wias-berlin.de/software/index.jsp?id=TetGen&lang=1)
- [pybind11](https://github.com/pybind/pybind11)
#### python
- [numpy](https://numpy.org)
#### build
- [cmake](https://cmake.org)
- [setuptools](https://setuptools.pypa.io/en/latest/)
- [wheel](https://wheel.readthedocs.io/en/stable/)
- [cibuildwheel](https://github.com/pypa/cibuildwheel)
