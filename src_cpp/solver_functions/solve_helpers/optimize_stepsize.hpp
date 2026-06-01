#pragma once
#include "../../common.hpp"

double optimize_stepsize(
    const Params& params,
    const Paths& paths,
    const SolutionData& data,
    const Field& gradient
);
