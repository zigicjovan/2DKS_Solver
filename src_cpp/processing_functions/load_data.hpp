#pragma once
#include "../common.hpp"

SolutionData load_solution_data(
    const Paths& paths,
    const Params& params
);

Field load_solution(const Paths& paths);
Field load_optimal_state(const Paths& paths);
std::vector<double> load_energy(const Paths& paths);
Diagnostics load_diagnostics(const Paths& paths);
