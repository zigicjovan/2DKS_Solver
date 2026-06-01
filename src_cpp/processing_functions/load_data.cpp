#include "load_data.h"
#include <fstream>

SolutionData load_solution_data(
    const Paths& paths,
    const Params& params
) {
    SolutionData data;

    data.u_hat = load_solution(paths);
    data.optimal_u_hat = load_optimal_state(paths);
    data.diag.energy = load_energy(paths);
    data.diag = load_diagnostics(paths);

    return data;
}

Field load_solution(const Paths& paths) {
    Field u;

    // TODO:
    // Replace with HDF5/binary loading.
    // For now, placeholder only.

    std::cout << "Loading solution from " << paths.solution_file << "\n";

    return u;
}

Field load_optimal_state(const Paths& paths) {
    Field u;

    std::cout << "Loading optimal state from " << paths.optimal_file << "\n";

    return u;
}

std::vector<double> load_energy(const Paths& paths) {
    std::vector<double> energy;

    std::ifstream file(paths.energy_file);
    double value;

    while (file >> value) {
        energy.push_back(value);
    }

    return energy;
}

Diagnostics load_diagnostics(const Paths& paths) {
    Diagnostics diag;

    // TODO:
    // Parse diagnostics file:
    // iteration, objective, gradient norm, step size, etc.

    return diag;
}
