#pragma once
#include "../../common.h"

double optimize_stepsize(
    const Params& params,
    const Paths& paths,
    const SolutionData& data,
    const Field& gradient
);
