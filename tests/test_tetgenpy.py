import unittest

import numpy as np

import tetgenpy


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


def points_and_tets(tetout):
    """
    Return points and tets from TetgenIO object
    """
    p = tetout.points()
    t = tetout.tetrahedra()
    return p, t


def has_points_and_tets(tetout):
    assert len(tetout.points())
    assert len(tetout.tetrahedra())


def same_or_more_points(original_points, tetout):
    assert len(tetout.points()) >= len(original_points)


def no_points_in_this_bound(points, bound):
    for i in range(len(bound[0])):
        in_range = np.logical_and(
            points[:, i] > bound[0][i], points[:, i] < bound[1][i]
        )
        assert sum(in_range) == 0, f"points detected in ({i})-dim"


def points_in_this_bound(points, bound):
    for i in range(len(bound[0])):
        in_range = np.logical_and(
            points[:, i] > bound[0][i], points[:, i] < bound[1][i]
        )
        assert sum(in_range) > 0, f"points missing in ({i})-dim"


def boxplc_using_holes():
    plc = tetgenpy.PLC()
    plc.add_points(vertices())
    plc.add_points(inner_vertices())
    # front
    # give none as hole.
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

    return plc


class tetgenpyTest(unittest.TestCase):
    """
    test to see if tets are generated with given input.
    """

    def test_plc2box_id_based(self):
        plc = tetgenpy.PLC()
        plc.add_points(vertices())
        plc.add_facets(faces())

        # call tetrahedralize without any switches
        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

    def test_plc2box_coord_based(self):
        plc = tetgenpy.PLC()
        plc.add_facets(vertices()[faces()])
        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

    def test_plc2box_with_a_hole_using_facet_with_holes(self):
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
            [
                [9, 11, 15, 13],
                [13, 15, 14, 12],
                [12, 14, 10, 8],
                [8, 10, 11, 9],
            ],
            4,
        )
        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

        tol = 1e-8
        iv = inner_vertices()
        lower_b = iv.min(axis=0)
        upper_b = iv.max(axis=0)

        no_points_in_this_bound(tetgenout.points(), [lower_b, upper_b])

    def test_plc2box_with_a_hole_using_holes(self):
        plc = boxplc_using_holes()

        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

        tol = 1e-8
        iv = inner_vertices()
        lower_b = iv.min(axis=0) + tol
        upper_b = iv.max(axis=0) - tol

        no_points_in_this_bound(tetgenout.points(), [lower_b, upper_b])

    def test_pl2longbox_with_hollow_half(self):
        plc = boxplc_using_holes()

        # attach another box along y-axis
        plc.add_points(
            vertices()[[2, 3, 6, 7]] + [0, 1, 0]
        )  # [16, 17, 18, 19]
        plc.add_facets(
            [
                [16, 17, 19, 18],
                [7, 19, 18, 6],
                [6, 18, 16, 2],
                [2, 16, 17, 3],
                [3, 17, 19, 7],
            ]
        )
        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

        tol = 1e-8
        iv = inner_vertices()
        lower_b = iv.min(axis=0) + tol
        upper_b = iv.max(axis=0) - tol

        no_points_in_this_bound(tetgenout.points(), [lower_b, upper_b])

        # extended box
        points_in_this_bound(
            tetgenout.points(),
            [[0 - tol, 1 - tol, 0 - tol], [1 + tol, 2 + tol, 1 + tol]],
        )

    def test_plc2long_box_with_2regions(self):
        plc = tetgenpy.PLC()
        plc.add_points(vertices())
        plc.add_facets(faces())
        plc.add_points(
            vertices()[[2, 3, 6, 7]] + [0, 1, 0]
        )  # [16, 17, 18, 19]
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
                [0.1, 0.1, 0.1, 10, 0.005],
                [0.1, 1.1, 0.1, 20, 1],
            ]
        )

        tetgenout = tetgenpy.tetrahedralize.tetrahedralize(
            "Aaq", plc.to_tetgenio()
        )

        has_points_and_tets(tetgenout)

        p = tetgenout.points()

        # overlap points in on y=1
        tol = 1e-8
        first_box_points = p[p[:, 1] < 1 + tol]
        second_box_points = p[p[:, 1] > 1 - tol]

        assert len(first_box_points) > len(
            second_box_points
        ), "add_regions() didn't work."


if __name__ == "__main__":
    unittest.main()
