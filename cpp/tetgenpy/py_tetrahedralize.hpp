#pragma once

#include <tetgen.h>

#include <pybind11/pybind11.h>

#include <tetgenpy/py_tetgenio.hpp>

namespace tetgenpy {

namespace py = pybind11;

void Tetrahedralize(std::string switches, PyTetgenIo& in, PyTetgenIo& out) {
  char* c_switches = switches.data();
  tetrahedralize(c_switches, &in, &out);
  in.is_output_ = false;
  out.is_output_ = true;
}

inline void add_pytetrahedralize(py::module_& m) {

  m.def("tetrahedralize", &Tetrahedralize);
}

} // namespace tetgenpy
