set(TETGENPY_SRCS
        init/py_tetgenio.cpp
        init/py_tetrahedralize.cpp
        tetgenpy_core.cpp)

# add target
pybind11_add_module(tetgenpy_core ${TETGENPY_SRCS})
target_include_directories(tetgenpy_core PRIVATE "${PROJECT_SOURCE_DIR}/cpp")
target_compile_definitions(tetgenpy_core PRIVATE ${TETGENPY_DEFS})
target_compile_options(tetgenpy_core
    PRIVATE
        ${TETGENPY_FLAGS}
        ${TETGENPY_OPTIMIZATION_FLAGS})
target_compile_features(tetgenpy_core PUBLIC cxx_std_11)

# link splinepy
target_link_libraries(tetgenpy_core PRIVATE tetgen)

# configure config files
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${cfg_dest}"
)

install(
    DIRECTORY "${PROJECT_SOURCE_DIR}/cpp/tetgenpy"
    DESTINATION ${incl_dest}
    FILES_MATCHING PATTERN "*.hpp"
)
