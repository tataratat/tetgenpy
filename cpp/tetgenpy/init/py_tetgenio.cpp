#include <tetgenpy/py_tetgenio.hpp>

namespace tetgenpy {
namespace py = pybind11;
void init_pytetgenio(py::module_& m) {add_pytetgenio_class(m);}
}
