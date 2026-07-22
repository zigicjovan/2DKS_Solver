#ifndef SOLVER_H
#define SOLVER_H

#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "FFTWPlanner.h"
#include "Timer.h"
#include "MPIContext.h"

#include <complex>

using namespace std;

enum StateSolutionType {
    SolveInitialState,               
    SolveTerminalState,    
    SolveBackwardInitialState,  
};

enum TimeSolutionType {
    SolveForwardInTime,       
    SolveBackwardInTime,        
};

enum OptimizeSolutionType { 
    OptimizeEnergyAmplification,      
    OptimizeLineSearchStepSize,           
};

struct EnergyData {
    double dEnergy;
    double dTimepoint;
};

class Solver{
private:
    Parameters& _params;
    Pathnames& _paths;
    FFTWPlanner& _fftwPlan;
    Timer& _timer;
    const MPIContext& _mpi;

    void saveSolutionDiagnostics(const vector<array<double, 4>>& vDiagnostics);
    void saveSolutionSpectrum(const vector<vector<double>>& vSpectrumHistory);
    void saveOptimizationDiagnostics(const array<double, 7>& vDiagnostics);
    void saveSolutionBranchAndPowerLaws(const array<double, 7>& vOptimalEnergySolution);
    void saveLineSearch(vector<double>& vLineSearchHistory);
    EnergyData getMaxEnergyL2InTimeWindow();

    void checkCFL(SolutionData& vData1, SolutionData& vData2);
    void checkSpectralResolution(const SolutionData& state);
    void setInitialCondition(SolutionData& vTargetState);
    void findContinuationForInitialData(SolutionData& vTargetState);
    void solveForwardInTime(SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
    void saveForwardState(double dTimePoint, size_t forwardIndex, size_t fullSteps, size_t stepsPerFile,  size_t totalSteps,
                          SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vStateCurrent);
    void loadForwardState(double dTimePoint, size_t forwardIndex, size_t fullSteps, size_t stepsPerFile, 
                          SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vForwardStateCurrent);
    void solveBackwardInTime(SolutionData& vObjectiveGradient, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, const SolutionData& vTargetEnd);
    void solveRiemmanianOptimization(double& dTargetValue, SolutionData& vObjectiveGradient, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                     SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
    void solveLineSearchOptimization(double& dTargetValue, SolutionData& DirectionCurr, SolutionData& vTargetStart, 
                                     SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
public:
    Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, const MPIContext& mpi);

    void setSolutionState(StateSolutionType targetType, SolutionData& vTargetState);
    void setSolutionInTime(TimeSolutionType targetType, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                           SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
    double getOptimalSolution(OptimizeSolutionType targetType, SolutionData& vObjectiveGradient, 
                              SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd);
};

#endif