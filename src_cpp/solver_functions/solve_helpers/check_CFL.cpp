#include "check_CFL.hpp"

void check_CFL(const Params& params) {
    double dx = params.Lx / params.Nx;
    double dy = params.Ly / params.Ny;

    double h = std::min(dx, dy);

    // Placeholder CFL-like diagnostic.
    //
    // For IMEX KS, the stiff linear part is treated implicitly,
    // so this is not the usual explicit heat-equation restriction.
    // Still useful as a sanity warning for nonlinear/advection-like terms.

    double cfl_indicator = params.dt / h;

    if (cfl_indicator > 1.0) {
        std::cerr
            << "Warning: CFL indicator dt/min(dx,dy) = "
            << cfl_indicator
            << " may be large.\n";
    }
}
