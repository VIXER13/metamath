cmake_minimum_required(VERSION 3.10)
project(metamath_example)

add_subdirectory(metamath)

set(CMAKE_CXX_STANDARD 17)

add_compile_options(-O2)

add_executable(metamath_example main.cpp)
target_link_libraries(metamath_example metamath_lib)
#target_link_libraries(metamath_example metamath_lib mesh_lib finite_element_solvers_lib)