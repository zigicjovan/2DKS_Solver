#include "solve_PDE.hpp"

#include "solve_helpers/initialcondition.hpp"
#include "solve_helpers/check_CFL.hpp"
#include "../processing_functions/save_data.hpp"

SolutionData solve_PDE(
    const Params& params,
    const Paths& paths,
    SolutionData data
) {
    const int N = params.Nx * params.Ny;

    // Allocate if empty.
    if (data.u_hat.empty()) {
        data.u_hat.resize(N);
        data.u_hat = initialcondition(params);
    }

    // Domain/operator setup.
    // TODO:
    // - build kx, ky arrays
    // - build k2, k4
    // - build dealias mask
    // - build IMEX RK coefficients
    //
    // This is where your MATLAB spectral operator setup should go.

    check_CFL(params);

    if (params.solve_forward) {
        std::cout << "Solving forward PDE...\n";

        for (int n = 0; n < params.nsteps; ++n) {
            // TODO:
            // Replace this with your MEX forward time-stepper logic.
            //
            // Expected future form:
            //
            // ks_forward_step(
            //     data.u_hat,
            //     operators,
            //     workspace,
            //     params.dt
            // );
            //
            // Important:
            // - FFTW inverse transforms require explicit division by Nx*Ny.
            // - Make sure wavenumber ordering matches MATLAB.
            // - Preserve MATLAB/MEX indexing until validation is complete.

            if (n % 100 == 0) {
                double energy = 0.0;

                for (const auto& z : data.u_hat) {
                    energy += std::norm(z);
                }

                data.diag.energy.push_back(energy);
            }
        }

        if (params.save_full_solution) {
            save_solution(paths, data.u_hat);
        }
    }

    if (params.solve_backward) {
        std::cout << "Solving backward/adjoint PDE...\n";

        // TODO:
        // Insert adjoint/backward solve here.
        // This should probably share the same operator setup as forward solve.
    }

    return data;
}
