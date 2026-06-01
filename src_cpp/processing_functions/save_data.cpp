#include "save_data.hpp"
#include <fstream>

void save_solution(const Paths& paths, const Field& u_hat) {
    std::ofstream file(paths.solution_file, std::ios::binary);

    std::size_t n = u_hat.size();
    file.write(reinterpret_cast<const char*>(&n), sizeof(std::size_t));
    file.write(
        reinterpret_cast<const char*>(u_hat.data()),
        static_cast<std::streamsize>(n * sizeof(Complex))
    );
}

void save_optimal_state(const Paths& paths, const Field& u_hat) {
    std::ofstream file(paths.optimal_file, std::ios::binary);

    std::size_t n = u_hat.size();
    file.write(reinterpret_cast<const char*>(&n), sizeof(std::size_t));
    file.write(
        reinterpret_cast<const char*>(u_hat.data()),
        static_cast<std::streamsize>(n * sizeof(Complex))
    );
}

void save_energy(const Paths& paths, const std::vector<double>& energy) {
    std::ofstream file(paths.energy_file);

    for (double e : energy) {
        file << e << "\n";
    }
}

void save_spectrum(const Paths& paths, const std::vector<double>& spectrum) {
    std::ofstream file(paths.spectrum_file);

    for (double s : spectrum) {
        file << s << "\n";
    }
}

void save_diagnostics(const Paths& paths, const Diagnostics& diag) {
    std::ofstream file(paths.diagnostics_file);

    file << "# iter gradient_norm step_size objective\n";

    std::size_t n = diag.gradient_norm.size();

    for (std::size_t i = 0; i < n; ++i) {
        double g = diag.gradient_norm[i];

        double a = i < diag.step_size.size()
            ? diag.step_size[i]
            : 0.0;

        double obj = i < diag.objective.size()
            ? diag.objective[i]
            : 0.0;

        file << i << " " << g << " " << a << " " << obj << "\n";
    }
}
