cmake_minimum_required(VERSION 3.12.0)

project(laplace_benchmark VERSION 1.0)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED YES)

add_compile_options(
  "-Wall" "-Wpedantic" "-Wextra" "-fexceptions"
  "$<$<CONFIG:DEBUG>:-O0;-g>"
)

add_library(laplace_benchmark SHARED src/laplace_2d.cpp)
set_target_properties(laplace_benchmark PROPERTIES PUBLIC_HEADER include/laplace_benchmark.h)
set_target_properties(laplace_benchmark PROPERTIES VERSION ${PROJECT_VERSION})
target_include_directories(laplace_benchmark PRIVATE include)
target_include_directories(laplace_benchmark PRIVATE src)

find_package(OpenMP REQUIRED)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")

option(COMPILE_PYBIND_MODULES "If turn ON will compile the library interface for python3" OFF)
if(COMPILE_PYBIND_MODULES)

  find_package(pybind11 REQUIRED)
  pybind11_add_module(pylbench python/py_laplace_benchmark.cpp)
  target_link_libraries(pylbench PUBLIC pybind11::module)
  target_include_directories(pylbench PRIVATE include)
  target_link_libraries(pylbench PRIVATE laplace_benchmark)

endif()

option(COMPILE_EXAMPLES "Set whether the examples are compiled" OFF)
if(COMPILE_EXAMPLES)

    set(EXAMPLES_INSTALL_PATH "${CMAKE_INSTALL_PREFIX}/bin" CACHE PATH "Set example installation path")

    add_executable(multithread_laplace_benchmark examples/multithread_benchmark.cpp)
    target_include_directories(multithread_laplace_benchmark PRIVATE include)
    target_link_libraries(multithread_laplace_benchmark PRIVATE laplace_benchmark)

    install(
        TARGETS multithread_laplace_benchmark
        RUNTIME DESTINATION "${EXAMPLES_INSTALL_PATH}"
        )
    
endif()

install(
    TARGETS laplace_benchmark
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_PREFIX}/include"
    )