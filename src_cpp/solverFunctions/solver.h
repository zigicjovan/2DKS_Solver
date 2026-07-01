#ifndef SOLVER_H
#define SOLVER_H

#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "FFTWPlanner.h"
#include "Timer.h"

#include <complex>

enum SolutionType {
    // TO DO: setSolutionState: initialize/load/save optIC, save {bwd,optTC (every time)}
    SolveInitialState,               
    SolveTerminalState,  
    // TO DO: setSolutionInTime: load data, check CFL, Nonlinear dealiasing, IMEX RK4, save {fwdALL (optimizer)}, save {energy,spectrum (last time)}, load {fwd (bwd)}  
    SolveForwardInTime,       
    SolveBackwardInTime, 
    // TO DO: optimizeSolution: RCG, brent, save {optIC, branch, diags, brent} 
    OptimizeEnergyAmplification,      
    OptimizeLineSearchStepSize,           
};

using Complex = std::complex<double>;

void setInitialCondition(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, SolutionData& vTargetState);
void setSolutionState(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, SolutionType targetType, SolutionData& vTargetState);
void setSolutionInTime(const Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, SolutionType targetType, SolutionData& vTargetStart, 
    SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
double getOptimalSolution(const Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, SolutionType targetType, SolutionData& vObjectiveGradient, 
    SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);

#endif