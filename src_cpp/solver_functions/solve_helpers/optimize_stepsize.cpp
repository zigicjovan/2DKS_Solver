#include "optimize_stepsize.hpp"
#include "../solve_PDE.hpp"

static double objective_function(
    const Params& params,
    const Paths& paths,
    const SolutionData& base_data,
    const Field& gradient,
    double alpha
) {
    SolutionData trial = base_data;

    for (std::size_t i = 0; i < trial.u_hat.size(); ++i) {
        trial.u_hat[i] += alpha * gradient[i];
    }

    Params forward_params = params;
    forward_params.solve_forward = true;
    forward_params.solve_backward = false;
    forward_params.save_full_solution = false;

    trial = solve_PDE(forward_params, paths, trial);

    if (trial.diag.energy.empty()) {
        return 0.0;
    }

    // Maximize terminal energy.
    return trial.diag.energy.back();
}

double optimize_stepsize(
    const Params& params,
    const Paths& paths,
    const SolutionData& data,
    const Field& gradient
) {
    // First-draft Brent placeholder.
    //
    // TODO:
    // Replace with full Brent maximization.
    // For now, test a few trial alphas and take the best.

    std::vector<double> candidates = {
        0.0,
        1e-4,
        1e-3,
        1e-2,
        1e-1
    };

    double best_alpha = 0.0;
    double best_obj = -1.0e300;

    for (double alpha : candidates) {
        double obj = objective_function(
            params,
            paths,
            data,
            gradient,
            alpha
        );

        if (obj > best_obj) {
            best_obj = obj;
            best_alpha = alpha;
        }
    }

    return best_alpha;
}
