#pragma once

#include <complex>
#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include <chrono>

using Complex = std::complex<double>;
using Field   = std::vector<Complex>;

enum class RunMode {
    ForwardOnly,
    Optimize
};

struct Paths {
    std::string output_dir      = "output";
    std::string solution_file   = "output/solution.bin";
    std::string energy_file     = "output/energy.txt";
    std::string diagnostics_file = "output/diagnostics.txt";
    std::string spectrum_file   = "output/spectrum.txt";
    std::string optimal_file    = "output/optimal_state.bin";
};

struct Params {
    // Grid/domain
    int Nx = 128;
    int Ny = 128;
    double Lx = 2.0 * M_PI;
    double Ly = 2.0 * M_PI;

    // Time stepping
    double dt = 1e-4;
    double T  = 1.0;
    int nsteps = 10000;

    // PDE parameters
    double ic_magnitude = 1.0;

    // Optimization
    int max_opt_iters = 20;
    double grad_tol = 1e-8;
    double brent_tol = 1e-6;

    // Runtime controls
    RunMode mode = RunMode::Optimize;
    bool continue_from_file = false;
    bool save_full_solution = true;
    bool delete_solution_after_processing = false;
    bool solve_forward = true;
    bool solve_backward = false;
};

struct Diagnostics {
    std::vector<double> energy;
    std::vector<double> gradient_norm;
    std::vector<double> step_size;
    std::vector<double> objective;
};

struct SolutionData {
    Field u_hat;          // Current spectral state
    Field optimal_u_hat;  // Best/optimized state
    Diagnostics diag;
};
