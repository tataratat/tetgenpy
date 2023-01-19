#pragma once

#include <algorithm>

#include <tetgen.h>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <tetgenpy/print.hpp>

namespace tetgenpy {

namespace py = pybind11;

template<typename ValueType>
static bool CheckPyArrayShape(const py::array_t<ValueType> arr,
                              const std::vector<int>& shape,
                              const bool throw_ = true) {
  const std::size_t expected_dim = shape.size();
  if (expected_dim != static_cast<std::size_t>(arr.ndim())) {
    if (!throw_)
      return false;
    PrintAndThrowError("Array dim mismatch.",
                       "Expected -",
                       expected_dim,
                       "Given -",
                       arr.ndim());
  }
  const py::ssize_t* arrshape = arr.shape();
  for (std::size_t i{}; i < expected_dim; ++i) {
    const int& shape_i = shape[i];
    if (shape_i < 0) {
      continue;
    } else {
      if (shape_i != arrshape[i]) {
        if (!throw_)
          return false;
        PrintAndThrowError("Array shape mismatch",
                           "in dimension [",
                           i,
                           "].",
                           "Expected -",
                           shape_i,
                           "Given -",
                           arrshape[i]);
      }
    }
  }
  return true;
}

class PyTetgenIo : public tetgenio {

public:
  using Base_ = tetgenio;

  // define constants
  static const int dim_{3};
  static const int n_region_entries_{5};
  static const int n_tetrahedron_vertices_order1_{4};
  static const int n_tetrahedron_vertices_order2_{10};

  // default ctor initializes. dtor clears.
  PyTetgenIo() { Base_::initialize(); }

  /// load full set of input
  /// split input facets into 2 types:
  ///   1) facets without a hole
  ///   2) facets with holes (we will refer to this as h_facet)
  ///
  /// points: (a, 3) np.ndarray
  ///   all points that appear for tetgen query
  /// facets: list
  ///   coplanar single polygons
  /// facet_markers: (len(facets),) np.ndarray
  ///   facet marker. If this facet is at boundary, this will be boundary id.
  ///   len(facet_markers) == len(facets) should hold
  /// h_facets: list
  ///   nested list for self contained facet description with holes list,
  ///   as well as marker.
  ///   [list of h_facet polygons, list of holes, facet marker]
  /// holes: (b, 3) np.ndarray
  ///   list of holes. this pokes holes.
  /// regions: (c, 5) np.ndarray
  ///   list of regions.
  ///   coordinate to mark the region, region attribute, max volume
  void Setup(py::array_t<REAL> points,
             py::list facets,
             py::array_t<int> facet_markers,
             py::list h_facets,
             py::array_t<REAL> holes,
             py::array_t<REAL> regions,
             bool debug = false) {

    PrintDebug(debug, "Starting PyTetgenIo::Load");

    // define some constants
    const int n_polygon_per_facet{1};

    /* numberofpoints, pointlist */
    // points shape test - 3D!
    CheckPyArrayShape(points, {-1, dim_});

    Base_::numberofpoints = static_cast<int>(points.shape(0));
    PrintDebug(debug, "set numberofpoints:", Base_::numberofpoints);

    Base_::pointlist = new REAL[Base_::numberofpoints * dim_];
    std::copy_n(static_cast<REAL*>(points.request().ptr),
                Base_::numberofpoints * dim_,
                Base_::pointlist);
    PrintDebug(debug, "set pointlist.");

    /* numberoffacets, facetlist, facetmarkerlist */
    const int n_facets = static_cast<int>(facets.size());
    const int n_facet_markers = static_cast<int>(facet_markers.size());
    const int n_h_facets = static_cast<int>(h_facets.size());
    Base_::numberoffacets = n_facets + n_h_facets;
    Base_::facetlist = new Base_::facet[Base_::numberoffacets];
    PrintDebug(debug, "set numberoffacets:", Base_::numberoffacets);
    PrintDebug(debug, "-> number of regular facets:", n_facets);
    PrintDebug(debug, "-> number of facets with holes:", n_h_facets);

    // following sizes should match
    // facet markers for h_facets will be added later
    if (n_facet_markers != 0 && n_facets != n_facet_markers) {
      PrintAndThrowError("facets and facet_markers should have same size.",
                         "facets:",
                         n_facets,
                         "/ facet_markers:",
                         n_facet_markers);
    }

    // loop and fillout facets
    // pointers for easy handling of raw ptr elements
    Base_::facet* f;
    Base_::polygon* p;
    for (int i{}; i < n_facets; ++i) {

      // TODO does this work without any issue? else, try reinterpret_steal
      const auto polygon_connec = facets[i].cast<py::list>();

      // get facet and init
      f = &Base_::facetlist[i];
      Base_::init(f);

      // fill polygon
      f->numberofpolygons = n_polygon_per_facet;
      f->polygonlist = new Base_::polygon[n_polygon_per_facet];
      p = &f->polygonlist[0]; // n_polygon_per_facet is 1
      Base_::init(p);
      const int n_polygon_vertices = static_cast<int>(polygon_connec.size());
      p->numberofvertices = n_polygon_vertices;
      p->vertexlist = new int[n_polygon_vertices];

      int j{};
      for (py::handle vertex_id : polygon_connec) {
        p->vertexlist[j] = vertex_id.cast<int>();
        ++j;
      }
    }
    PrintDebug(debug, "set facetlist.");

    // set facet markers. could put it to the loop above,
    // but maybe not everyone needs this
    if (n_facet_markers != 0) {
      Base_::facetmarkerlist = new int[Base_::numberoffacets];
      std::copy_n(static_cast<int*>(facet_markers.request().ptr),
                  // Base_::numberoffacets, //
                  n_facets, //
                  Base_::facetmarkerlist);
      PrintDebug(debug, "set facetmarkerlist for facets.");
    }

    // fill h_facets
    int h{n_facets}; // hfacet counter
    for (py::handle hf : h_facets) {
      // index hints
      // 0. -> list of list (polygons; int)
      // 1. -> list of list (holes; double)
      // 2. -> facet marker (;int)
      const auto hf_description = hf.cast<py::list>();
      const auto list_of_polygons = hf_description[0].cast<py::list>();
      const auto list_of_holes = hf_description[1].cast<py::list>();
      const int facet_mark = hf_description[2].cast<int>();

      f = &Base_::facetlist[h];
      Base_::init(f);

      const int n_poly_per_hf = static_cast<int>(list_of_polygons.size());
      f->numberofpolygons = n_poly_per_hf;
      f->polygonlist = new Base_::polygon[n_poly_per_hf];

      // 0
      int i{}; // polygon counter
      for (py::handle poly : list_of_polygons) {
        p = &f->polygonlist[i];
        Base_::init(p);
        const auto poly_vertexlist = poly.cast<py::list>();
        const int n_vertices_per_poly =
            static_cast<int>(poly_vertexlist.size());
        p->numberofvertices = n_vertices_per_poly;
        p->vertexlist = new int[n_vertices_per_poly];

        // fill vertices
        int j{};
        for (py::handle vertex_id : poly_vertexlist) {
          p->vertexlist[j] = vertex_id.cast<int>();
          ++j;
        }
        ++i;
      }

      // 1
      const int n_hf_holes = static_cast<int>(list_of_holes.size());
      f->numberofholes = n_hf_holes;
      f->holelist = new REAL[n_hf_holes * dim_];
      i = 0;
      for (py::handle hole : list_of_holes) {
        // this needs to have len 3
        int j{};
        const auto hole_coord = hole.cast<py::list>();
        if (hole_coord.size() != 3) {
          PrintAndThrowError("h_facet holes must be 3D coordinates.",
                             "Given hole has the size of (",
                             hole_coord.size(),
                             ")");
        }
        for (py::handle coord : hole) {
          f->holelist[i * dim_ + j] = coord.cast<REAL>();
          ++j;
        }
        ++i;
      }

      // 2
      Base_::facetmarkerlist[h] = facet_mark;

      // next h_facet
      ++h;
    }

    /* numberofholes, holelist */
    /* numberofregions, regionlist */
  }

  // output
  py::array_t<REAL> GetPoints() {
    if (pointlist != (REAL*) NULL) {
      py::array_t<REAL> points(Base_::numberofpoints * dim_);
      std::copy_n(Base_::pointlist,
                  Base_::numberofpoints * dim_,
                  static_cast<REAL*>(points.request().ptr));
      points.resize({Base_::numberofpoints, dim_});
      return points;
    }
    // no points
    return py::array_t<REAL>(0);
  }

  py::array_t<int> GetTetrahedrons() {
    if (tetrahedronlist != (int*) NULL) {
      // order 1 for now
      py::array_t<int> tets(Base_::numberoftetrahedra
                            * n_tetrahedron_vertices_order1_);
      std::copy_n(Base_::tetrahedronlist,
                  Base_::numberoftetrahedra * n_tetrahedron_vertices_order1_,
                  static_cast<int*>(tets.request().ptr));
      tets.resize({Base_::numberoftetrahedra, n_tetrahedron_vertices_order1_});
      return tets;
    }
    // no tets
    return py::array_t<int>(0);
  }
};

inline void add_pytetgenio_class(py::module_& m) {

  py::class_<PyTetgenIo> klasse(m, "tetgenio");

  klasse.def(py::init<>())
      .def("setup",
           &PyTetgenIo::Setup,
           py::arg("points"),
           py::arg("facets"),
           py::arg("facet_markers"),
           py::arg("h_facets"),
           py::arg("holes"),
           py::arg("regions"),
           py::arg("debug"))
      .def("points", &PyTetgenIo::GetPoints)
      .def("tets", &PyTetgenIo::GetTetrahedrons);
}

} // namespace tetgenpy
