#ifndef SOLVER_H
#define SOLVER_H

#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "FFTWPlanner.h"
#include "Timer.h"

#include <complex>

using namespace std;

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

class Solver{
private:
    Parameters& _params;
    Pathnames& _paths;
    FFTWPlanner& _fftwPlan;
    Timer& _timer;

    void saveSolutionDiagnostics(const vector<array<double, 4>>& vDiagnostics);
    void saveSolutionSpectrum(const vector<vector<double>>& vSpectrumHistory);
    void checkCFL(SolutionData& vData1, SolutionData& vData2);
    void setInitialCondition(SolutionData& vTargetState);
public:
    Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer);

    void setSolutionState(StateSolutionType targetType, SolutionData& vTargetState);
    void setSolutionInTime(TimeSolutionType targetType, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
        SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
    double getOptimalSolution(OptimizeSolutionType targetType, SolutionData& vObjectiveGradient, 
        SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
};

#endif