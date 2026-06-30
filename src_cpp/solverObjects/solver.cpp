#include "solver.h"
#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"

void setSolutionState(const Parameters& params, const Pathnames& paths, SolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        case SolveInitialState: {
            break;
        }
        case SolveTerminalState: {
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break;
    }  
}

void setSolutionInTime(const Parameters& params, const Pathnames& paths, SolutionType targetType, SolutionData& vTargetStart, 
    SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case SolveForwardInTime: {
            break;
        }
        case SolveBackwardInTime: {
            break;
        }
        case SolveInitialState:
        case SolveTerminalState:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break;
    }  
}

double getOptimalSolution(const Parameters& params, const Pathnames& paths, SolutionType targetType, SolutionData& vObjectiveGradient, 
    SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case OptimizeEnergyAmplification: {
            break;
        }
        case OptimizeLineSearchStepSize: {
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case SolveInitialState:
        case SolveTerminalState:
            break;
    } 
    return 0.0;
}
