#include "common.hpp"

#include "processing_functions/load_data.hpp"
#include "processing_functions/save_data.hpp"
#include "processing_functions/delete_data.hpp"

#include "solver_functions/solve_PDE.hpp"
#include "solver_functions/optimize_solution.hpp"

int main() {
    auto start_time = std::chrono::high_resolution_clock::now();

    Paths paths;
    Params params;

    // Main parameter settings.
    params.Nx = 128;
    params.Ny = 128;
    params.dt = 1e-4;
    params.T  = 1.0;
    params.nsteps = static_cast<int>(params.T / params.dt);

    params.mode = RunMode::Optimize;
    params.continue_from_file = false;

    SolutionData data;

    // Optional continuation from previous run.
    if (params.continue_from_file) {
        data = load_solution_data(paths, params);
    } else {
        // Initialize and solve forward once.
        data = solve_PDE(params, paths, data);
    }

    // Optimization loop.
    if (params.mode == RunMode::Optimize) {
        data = optimize_solution(params, paths, data);
    }

    // Save main scalar outputs.
    save_energy(paths, data.diag.energy);
    save_diagnostics(paths, data.diag);

    // Save optimized state if available.
    save_optimal_state(paths, data.optimal_u_hat);

    // Optional cleanup of heavy solution files.
    if (params.delete_solution_after_processing) {
        delete_solution(paths);
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    double elapsed =
        std::chrono::duration<double>(end_time - start_time).count();

    std::cout << "Run completed in " << elapsed << " seconds.\n";

    return 0;
}
