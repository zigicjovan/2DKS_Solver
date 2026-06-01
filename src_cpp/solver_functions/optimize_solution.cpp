#include "optimize_solution.h"

#include "solve_PDE.h"
#include "solve_helpers/optimize_stepsize.h"
#include "../processing_functions/save_data.h"

SolutionData optimize_solution(
    const Params& params,
    const Paths& paths,
    SolutionData data
) {
    std::cout << "Starting optimization...\n";

    data.optimal_u_hat = data.u_hat;

    for (int iter = 0; iter < params.max_opt_iters; ++iter) {
        std::cout << "Optimization iteration " << iter << "\n";

        // 1. Backward/adjoint solve.
        Params backward_params = params;
        backward_params.solve_forward = false;
        backward_params.solve_backward = true;

        data = solve_PDE(backward_params, paths, data);

        // 2. Compute gradient.
        Field gradient(data.u_hat.size());

        // TODO:
        // Replace this placeholder with adjoint-derived gradient.
        for (std::size_t i = 0; i < gradient.size(); ++i) {
            gradient[i] = data.u_hat[i];
        }

        double grad_norm = 0.0;
        for (const auto& g : gradient) {
            grad_norm += std::norm(g);
        }
        grad_norm = std::sqrt(grad_norm);

        data.diag.gradient_norm.push_back(grad_norm);

        if (grad_norm < params.grad_tol) {
            std::cout << "Gradient tolerance reached.\n";
            break;
        }

        // 3. Brent line search / step-size optimization.
        double alpha = optimize_stepsize(params, paths, data, gradient);
        data.diag.step_size.push_back(alpha);

        // 4. Apply gradient step.
        for (std::size_t i = 0; i < data.u_hat.size(); ++i) {
            data.u_hat[i] += alpha * gradient[i];
        }

        // 5. Forward solve with updated IC.
        Params forward_params = params;
        forward_params.solve_forward = true;
        forward_params.solve_backward = false;

        data = solve_PDE(forward_params, paths, data);

        // 6. Save diagnostics every iteration.
        save_diagnostics(paths, data.diag);
    }

    data.optimal_u_hat = data.u_hat;

    return data;
}
