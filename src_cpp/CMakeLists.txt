cmake_minimum_required(VERSION 3.16)

project(KS2D_CPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(ks2d
    src/main/main.cpp
    src/solver_functions/solve_PDE.cpp
    src/solver_functions/optimize_solution.cpp
    src/solver_functions/solve_helpers/initialcondition.cpp
    src/solver_functions/solve_helpers/optimize_stepsize.cpp
    src/solver_functions/solve_helpers/check_CFL.cpp
    src/processing_functions/load_data.cpp
    src/processing_functions/save_data.cpp
    src/processing_functions/delete_data.cpp
)

target_compile_options(ks2d PRIVATE -O3 -march=native -Wall -Wextra)

# Later add FFTW/OpenMP here:
# find_package(OpenMP REQUIRED)
# target_link_libraries(ks2d PRIVATE OpenMP::OpenMP_CXX fftw3 fftw3_threads m)
