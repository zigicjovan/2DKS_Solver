#include "delete_data.h"
#include <filesystem>

static void remove_if_exists(const std::string& filename) {
    if (std::filesystem::exists(filename)) {
        std::filesystem::remove(filename);
        std::cout << "Deleted " << filename << "\n";
    }
}

void delete_solution(const Paths& paths) {
    remove_if_exists(paths.solution_file);
}

void delete_optimal_state(const Paths& paths) {
    remove_if_exists(paths.optimal_file);
}

void delete_energy(const Paths& paths) {
    remove_if_exists(paths.energy_file);
}

void delete_spectrum(const Paths& paths) {
    remove_if_exists(paths.spectrum_file);
}

void delete_diagnostics(const Paths& paths) {
    remove_if_exists(paths.diagnostics_file);
}
