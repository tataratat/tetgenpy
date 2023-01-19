#include <pybind11/pybind11.h>

namespace tetgenpy {
namespace py = pybind11;
void init_pytetgenio(py::module_&);
} // namespace tetgenpy

PYBIND11_MODULE(tetgenpy_core, m) { tetgenpy::init_pytetgenio(m); }
