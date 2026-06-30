#include "readParameterSettings.h"
#include "writePathNames.h"
#include "Timer.h"
#include "solver.h"
#include "SolutionData.h"
#include "updateDirectory.h"

#include "iostream"

int main(int argc, char* argv[]) {
    // Step 1: prepare test case
    Parameters params = readParameterSettings(argc, argv);
    Pathnames paths = writePathNames(params);

    // Step 2: track computation time
    Timer timer;
    timer.start(); 
    // timer.printInterval("");

    // Step 3: allocate storage for states and time-dependent solutions
    SolutionData vStateInitial(params, InitialState); 
    SolutionData vStateTerminal(params, TerminalState); 
    SolutionData vObjectiveGradient(params, BackwardInitialState); 
    SolutionData vHistoryIntermediate(params, IntermediateHistory);
    SolutionData vHistoryRemainder(params, RemainderHistory);

    // Step 4: solve computational problem in Fourier domain using pseudo-spectral time-stepping method
    setSolutionState(params, paths, SolveInitialState, vStateInitial);
    setSolutionInTime(params, paths, SolveForwardInTime, vStateInitial, vHistoryIntermediate, vHistoryRemainder, vStateTerminal);
        // then if params.bOptimizeSolution == 1 updateDirectoryData
    if (params.bOptimizeSolution == 1) {
        double dOptimalSolution = getOptimalSolution(params, paths, OptimizeEnergyAmplification, vObjectiveGradient, vStateInitial, vHistoryIntermediate, vHistoryRemainder, vStateTerminal);  
            // then updateDirectoryData: save diags and delete forward solution in dir
        std::cout << dOptimalSolution;
    }

    // Step 5: Clean up  
    saveData(paths, vStateInitial, InitialState);
    // saveData(paths, vHistoryIntermediate, RemainderHistory, params.dTimeWindow);
    // loadData(paths, vHistoryRemainder, RemainderHistory);
    // deleteData(paths);
    timer.stop();

    return 0;
}