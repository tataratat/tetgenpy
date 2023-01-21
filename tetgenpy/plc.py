"""
Class to prepare input Piecewise Linear Complexes for tetgen.
"""


class PLC:
    def __init__(self):
        self._points = []
        self._facets = []
        self._facet_with_holes = []
        self._holes = []
        self._regions = []
        self._facet_markerts = []

        self._point_id_offset = 0

    def add_points(self, points):
        pass

    def add_facets(self, polygons, facet_id=-1):
        pass

    def add_facet_with_holes(self, polygons, holes, facet_id=-1):
        pass

    def add_holes(self, holes):
        pass

    def add_regions(self, regions):
        pass

    def sofarsogood(self):
        """sanity check"""
        pass

    def to_tetgenio(self):
        pass

    def show(self):
        pass
