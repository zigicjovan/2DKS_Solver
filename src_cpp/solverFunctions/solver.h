#ifndef SOLVER_H
#define SOLVER_H

#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "FFTWPlanner.h"
#include "Timer.h"

#include <complex>

enum StateSolutionType {
    // TO DO: setSolutionState: initialize/load/save optIC, save {bwd,optTC (every time)}
    SolveInitialState,               
    SolveTerminalState,      
};

enum TimeSolutionType {
    // TO DO: setSolutionInTime: load data, check CFL, Nonlinear dealiasing, IMEX RK4, save {fwdALL (optimizer)}, save {energy,spectrum (last time)}, load {fwd (bwd)}  
    SolveForwardInTime,       
    SolveBackwardInTime,        
};

enum OptimizeSolutionType {
    // TO DO: optimizeSolution: RCG, brent, save {optIC, branch, diags, brent} 
    OptimizeEnergyAmplification,      
    OptimizeLineSearchStepSize,           
};

using Complex = std::complex<double>;

void saveSolutionDiagnostics(const Parameters& params, const Pathnames& paths, const std::vector<std::array<double, 4>>& vDiagnostics);
void saveSolutionSpectrum(const Parameters& params, const Pathnames& paths, const std::vector<std::vector<double>>& vSpectrumHistory);
void checkCFL(const Parameters& params, SolutionData& vData1, SolutionData& vData2);
void setInitialCondition(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, SolutionData& vTargetState);
void setSolutionState(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, StateSolutionType targetType, SolutionData& vTargetState);
void setSolutionInTime(const Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, TimeSolutionType targetType, SolutionData& vTargetStart, 
    SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
double getOptimalSolution(Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, OptimizeSolutionType targetType, SolutionData& vObjectiveGradient, 
    SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);

#endif