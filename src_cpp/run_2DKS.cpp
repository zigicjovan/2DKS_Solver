#include "Parameters.h"
#include "Pathnames.h"
#include "FFTWPlanner.h"
#include "Timer.h"
#include "Solver.h"
#include "SolutionData.h"

#include <iostream>

int main(int argc, char* argv[]) {
    // Step 1: prepare test case
    Parameters params(argc, argv);
    Pathnames paths(params);
    FFTWPlanner fftwPlan(params);

    // Step 2: track computation time
    Timer timer;
    timer.start(); 

    // Step 3: allocate storage for states and time-dependent solutions
    SolutionData vStateInitial(params, InitialState); 
    SolutionData vStateTerminal(params, TerminalState); 
    SolutionData vObjectiveGradient(params, BackwardInitialState); 
    SolutionData vHistoryIntermediate(params, IntermediateHistory);
    SolutionData vHistoryRemainder(params, RemainderHistory);

    // Step 4: solve computational problem in Fourier domain using pseudo-spectral time-stepping method
    Solver solver(params, paths, fftwPlan, timer);
    solver.setSolutionState(SolveInitialState, vStateInitial);
    solver.setSolutionInTime(SolveForwardInTime, vStateInitial, vHistoryIntermediate, vHistoryRemainder, vStateTerminal);
    double dOptimalSolution = solver.getOptimalSolution(OptimizeEnergyAmplification, vObjectiveGradient, vStateInitial, 
                            vHistoryIntermediate, vHistoryRemainder, vStateTerminal);  
    std::cout << "Maximum Energy Amplification: " << dOptimalSolution;

    // Step 5: Clean up  
    // vStateInitial.deleteData(paths);
    timer.stop();

    return 0;
}