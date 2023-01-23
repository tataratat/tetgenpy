"""
Class to prepare input Piecewise Linear Complexes for tetgen.
"""
import numpy as np


def _check_2d_and_shape1(name, array, shape1):
    shape = array.shape
    if len(shape) != 2 or shape[1] != shape1:
        raise ValueError(
            f"{name} should be 2D array-like input with (n, {shape1}) shape."
        )


class PLC:
    def __init__(self):
        self._points = []
        self._facets = []
        self._facet_with_holes = []
        self._holes = []
        self._regions = []
        self._facet_markers = []

        # self._point_id_offset = 0
        self._facet_id_offset = 0

    def add_points(self, points):
        """
        Adds points. list, tuple, or np.ndarray.

        Parameters
        ----------
        points: (n, 3)  array-like

        Returns
        -------
        None
        """
        # keep points as np.ndarray
        points = np.asanyarray(points)

        _check_2d_and_shape1("points", points, 3)

        # append points
        self._points.append(points)

    def _process_polygons(self, polygons, coordinates=None):
        """
        Helper to process facets, in case they are given as coordinates form.
        Coordinates will be added to points using add_points() and they will
        be turned into int list that references added points.
        """
        if isinstance(polygons, np.ndarray):
            polygons = polygons.tolist()

        if not isinstance(polygons, (tuple, list)):
            raise TypeError(
                "facets/polygons should be tuple, list, or np.ndarray"
            )

        if coordinates is None:
            coordinates = isinstance(polygons[0][0], (tuple, list, np.ndarray))

        if coordinates:
            point_id_offset = len(self._points)
            id_based_polygons = []
            for i, coords in polygons:
                # add_points will check for valid input
                self.add_points(coords)

                # point to added points
                id_based_polygons.append(
                    list(range(point_id_offset, point_id_offset + len(coords)))
                )

            # at this point, we overwrite polygons
            polygons = id_based_polygons

        return polygons

    def add_facets(self, polygons, facet_id=-1, coordinates=None):
        """
        Adds facets with one polygon per facet.

        Parameters
        ----------
        polygons: array-like
          list of polygons by coordinate or indices referring to points list.
          coordinates are expected to have [[[x1_1, y1_2, z1_3], ...], ...]
          and indices [[p1_1, p1_2, ...], ... ].
        facet_id: int
          value used to set facetmarker.
          can be also interpreted as boundary id as all the output faces on
          this facet will have this id.
        coordinates: bool
          Default is None. Can specify if given polygons are coordinate or
          index based. If None, will check the first entry and determine.

        Returns
        -------
        None
        """
        # make sure polygons are list of list referring to points
        processed_polygon_ids = self._process_polygons(polygons, coordinates)

        # by now, we should have points and ids.
        self._facets.extend(processed_polygon_ids)
        self._facet_markers.extend([facet_id] * len(polygons))

    def add_facet_with_holes(
        self, polygons, holes, facet_id=-1, coordinates=None
    ):
        """
        Supports full features of tetgenio::facet.

        Parameters
        ----------
        polygons: array-like
          same as add_facets
        holes: array-like
          list of locations to poke holes in this facet. If you don't want
          holes, you can set it to None.
        facet_id: int
        coordinates: bool

        Returns
        -------
        None
        """
        # facets
        processed_polygon_ids = self._process_polygons(polygons, coordinates)

        facet_with_holes = [processed_polygon_ids]

        # holes
        if holes is not None:
            holes = np.asanyarray(holes)
            _check_2d_and_shape1("facet holes", holes)
            facet_with_holes.append(holes)
        else:
            # empty hole
            facet_with_holes.append([])

        # facet_id
        facet_with_holes.append(int(facet_id))

        self._facet_with_holes.append(facet_with_holes)

    def add_holes(self, holes):
        pass

    def add_regions(self, regions):
        pass

    def sofarsogood(self):
        """assert for sanity check"""

        assert isinstance(self._points, np.ndarray)
        assert self.points.ndim == 2
        assert self.points.shape[1] == 3

    def to_tetgenio(self):
        pass

    def show(self):
        pass
