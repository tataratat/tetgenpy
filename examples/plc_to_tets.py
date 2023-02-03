import numpy as np

import tetgenpy

has_gus = False
try:
    import gustaf as gus

    has_gus = True
except ImportError:
    pass


def vertices():
    """
    Unit cube vertices
    """
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ]
    )


def inner_vertices():
    """
    Unit cube hole  veritces
    """
    v = vertices()
    v *= [0.4, 1, 0.4]
    v += [0.3, 0, 0.3]
    return v


def faces():
    """
    Unit cube faces
    """
    return np.array(
        [
            [1, 0, 2, 3],
            [0, 1, 5, 4],
            [2, 0, 4, 6],
            [1, 3, 7, 5],
            [3, 2, 6, 7],
            [4, 5, 7, 6],
        ],
        dtype=np.int32,
    )


def gusshow(tetout):
    """
    shows output if `gustaf` is installed
    """
    if not has_gus:
        return None

    p = tetout.points()
    t = tetout.tetrahedra()
    gus.Volumes(p, t).shrink().show()


if __name__ == "__main__":
    # tetgen uses `tetgenio` objects to create input and output.
    # tetgenpy.core.TetgenIO is a derived class that can be directly used
    # for tetrahedralize().
    # or, you can also use a helper class PLC.
    plc = tetgenpy.PLC()

    # 1. points and facets (integers referring to id of the points)
    plc.add_points(vertices())
    plc.add_facets(faces())

    # call tetrahedralize without any switches
    tetgenout = tetgenpy.tetrahedralize("", plc.to_tetgenio())
    gusshow(tetgenout)

    # 2. define facets with sequence of coordiantes instead of point ids.
    plc = tetgenpy.PLC()
    plc.add_facets(vertices()[faces()])
    tetgenout = tetgenpy.tetrahedralize(
        "qa0.00005", plc.to_tetgenio()
    )
    gusshow(tetgenout)

    # 3. box with a hole using add_facet_with_holes
    plc = tetgenpy.PLC()
    plc.add_points(vertices())
    plc.add_points(inner_vertices())

    # front
    plc.add_facet_with_holes(
        polygons=[[0, 1, 5, 4], [8, 9, 13, 12]],
        holes=[[0.5, 0.0, 0.5]],
        facet_id=1,
    )
    # back
    plc.add_facet_with_holes(
        polygons=[[2, 3, 7, 6], [10, 11, 15, 14]],
        holes=[[0.5, 1.0, 0.5]],
        facet_id=2,
    )
    # outer sides
    plc.add_facets(
        [[1, 3, 7, 5], [5, 7, 6, 4], [4, 6, 2, 0], [0, 2, 3, 1]],
        3,
    )
    # inner sides
    plc.add_facets(
        [[9, 11, 15, 13], [13, 15, 14, 12], [12, 14, 10, 8], [8, 10, 11, 9]],
        4,
    )
    tetgenout = tetgenpy.tetrahedralize(
        "qa0.00005", plc.to_tetgenio()
    )
    gusshow(tetgenout)

    # 4. long box with "half" hole
    # this time, let's try it with add_holes()
    # for that, we can still use add_facet_with_holes(); just need to say
    # that there are no holes.
    # long box is create by attaching extra box.
    # that way we can also use add_regions() to assign different
    # volume constraints
    plc = tetgenpy.PLC()
    plc.add_points(vertices())
    plc.add_points(inner_vertices())
    # front
    # give None as hole.
    plc.add_facet_with_holes(
        [
            [0, 1, 5, 4],
            [8, 9, 13, 12],
        ],
        holes=None,
    )
    plc.add_facet_with_holes(
        [
            [2, 3, 7, 6],
            [10, 11, 15, 14],
        ],
        holes=None,
    )
    plc.add_facets(
        [
            [1, 3, 7, 5],
            [5, 7, 6, 4],
            [4, 6, 2, 0],
            [0, 2, 3, 1],
            [9, 11, 15, 13],
            [13, 15, 14, 12],
            [12, 14, 10, 8],
            [8, 10, 11, 9],
        ],
    )
    # this pokes holes at defined point and it propagates all its way to
    # the surrounding facets
    plc.add_holes([[0.5, 0.5, 0.5]])

    # attach another box
    plc.add_points(vertices()[[2, 3, 6, 7]] + [0, 1, 0])  # [16, 17, 18, 19]
    plc.add_facets(
        [
            [16, 17, 19, 18],
            [7, 19, 18, 6],
            [6, 18, 16, 2],
            [2, 16, 17, 3],
            [3, 17, 19, 7],
        ]
    )
    # [x,y,z,region attribute, volume constraint]
    plc.add_regions(
        [
            [0.1, 0.1, 0.1, 10, 0.00005],
            [0.1, 1.1, 0.1, 20, 0.1],
        ]
    )
    tetgenout = tetgenpy.tetrahedralize(
        "Aaq1.1/0", plc.to_tetgenio()
    )
    gusshow(tetgenout)
