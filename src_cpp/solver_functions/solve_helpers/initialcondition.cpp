#include "initialcondition.h"

Field initialcondition(const Params& params) {
    const int N = params.Nx * params.Ny;
    Field u_hat(N);

    // Placeholder spectral initial condition.
    //
    // TODO:
    // Replace with your MATLAB IC options:
    // - random noise
    // - random Fourier phases
    // - structured ICs: s, s1, sc1, mmhf, etc.
    //
    // Keep the same normalization convention as MATLAB.

    for (int j = 0; j < params.Ny; ++j) {
        for (int i = 0; i < params.Nx; ++i) {
            int idx = j * params.Nx + i;

            double value = 0.0;

            if (i == 1 && j == 0) {
                value = params.ic_magnitude;
            }

            u_hat[idx] = Complex(value, 0.0);
        }
    }

    return u_hat;
}
