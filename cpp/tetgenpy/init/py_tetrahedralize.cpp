#include <tetgenpy/py_tetrahedralize.hpp>

namespace tetgenpy {
namespace py = pybind11;
void init_pytetrahedralize(py::module_& m) { add_pytetrahedralize(m); }
} // namespace tetgenpy
