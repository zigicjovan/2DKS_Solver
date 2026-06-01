#pragma once
#include "../common.hpp"

void save_solution(const Paths& paths, const Field& u_hat);
void save_optimal_state(const Paths& paths, const Field& u_hat);
void save_energy(const Paths& paths, const std::vector<double>& energy);
void save_spectrum(const Paths& paths, const std::vector<double>& spectrum);
void save_diagnostics(const Paths& paths, const Diagnostics& diag);
