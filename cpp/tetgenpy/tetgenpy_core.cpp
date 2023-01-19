#include <pybind11/pybind11.h>

namespace tetgenpy {
namespace py = pybind11;
void init_pytetgenio(py::module_&);
void init_pytetrahedralize(py::module_&);
} // namespace tetgenpy

PYBIND11_MODULE(tetgenpy_core, m) {
  tetgenpy::init_pytetgenio(m);
  tetgenpy::init_pytetrahedralize(m);
}
