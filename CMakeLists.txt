cmake_minimum_required(VERSION 3.16)

project(KS2D_CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ks2d
    src_cpp/main/main.cpp
    src_cpp/solver_functions/solve_PDE/solve_PDE.cpp
    src_cpp/solver_functions/optimize_solution/optimize_solution.cpp
    src_cpp/solver_functions/solve_helpers/initialcondition/initialcondition.cpp
    src_cpp/solver_functions/solve_helpers/optimize_stepsize/optimize_stepsize.cpp
    src_cpp/solver_functions/solve_helpers/check_CFL/check_CFL.cpp
    src_cpp/processing_functions/load_data/load_data.cpp
    src_cpp/processing_functions/save_data/save_data.cpp
    src_cpp/processing_functions/delete_data/delete_data.cpp
)

target_compile_options(ks2d PRIVATE -O3 -march=native -Wall -Wextra)

# Later add FFTW/OpenMP here:
# find_package(OpenMP REQUIRED)
# target_link_libraries(ks2d PRIVATE OpenMP::OpenMP_CXX fftw3 fftw3_threads m)
