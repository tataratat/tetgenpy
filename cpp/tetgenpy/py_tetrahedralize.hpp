#pragma once

#include <tetgen.h>

#include <pybind11/pybind11.h>

#include <tetgenpy/py_tetgenio.hpp>

namespace tetgenpy {

namespace py = pybind11;

void Tetrahedralize(std::string switches, PyTetgenIo& in, PyTetgenIo& out) {
  char* c_switches = &switches.front();
  tetrahedralize(c_switches, &in, &out);
  in.is_output_ = false;
  out.is_output_ = true;
}

int Main(py::list sys_argv) {
  // prepare argc
  int argc = static_cast<int>(sys_argv.size());

  // prepare argv
  // TODO now sure how casting char* would work, so start from string
  std::vector<std::string> argv_string;
  argv_string.reserve(argc);
  for (auto& s_a : sys_argv) {
    argv_string.emplace_back(s_a.cast<std::string>());
  }

  std::vector<char*> argv_char;
  argv_char.reserve(argc);
  for (auto& a_s : argv_string) {
    argv_char.emplace_back(&a_s.front());
  }

  char** argv = argv_char.data();

  // following parts are direct copies from tetgen.cxx
  tetgenbehavior b;

  tetgenio in, addin, bgmin;

  if (!b.parse_commandline(argc, argv)) {
    terminatetetgen(NULL, 10);
  }

  // Read input files.
  if (b.refine) { // -r
    if (!in.load_tetmesh(b.infilename, (int) b.object)) {
      terminatetetgen(NULL, 10);
    }
  } else { // -p
    if (!in.load_plc(b.infilename, (int) b.object)) {
      terminatetetgen(NULL, 10);
    }
  }
  if (b.insertaddpoints) { // -i
    // Try to read a .a.node file.
    addin.load_node(b.addinfilename);
  }
  if (b.metric) { // -m
    // Try to read a background mesh in files .b.node, .b.ele.
    bgmin.load_tetmesh(b.bgmeshfilename, (int) b.object);
  }

  tetrahedralize(&b, &in, NULL, &addin, &bgmin);

  return 0;
}

inline void add_pytetrahedralize(py::module_& m) {

  m.def("tetrahedralize", &Tetrahedralize);
  m.def("main", &Main);
}

} // namespace tetgenpy
