#include "readParameterSettings.h"
#include "writePathNames.h"
#include "Timer.h"
#include "SolutionData.h"
// #include "updateDirectory.h"
// #include "computeSolution.h"
// #include "optimizeSolution.h"

int main(int argc, char* argv[]) {
    // Step 1: prepare test case
    Parameters params = readParameterSettings(argc, argv);
    Pathnames paths = writePathNames(params);

    // Step 2: track computation time
    Timer timer;
    timer.start(); 
    // timer.printInterval("");

    // Step 3: allocate storage for initial state and time-dependent solution
    SolutionData vStateInitial(params, SolutionDataType::InitialState); 
    SolutionData vHistoryIntermediate(params, SolutionDataType::IntermediateHistory);
    SolutionData vHistoryRemainder(params, SolutionDataType::RemainderHistory);
    // saveData(paths, vHistoryIntermediate, SolutionDataType::RemainderHistory, params.dTimeWindow);
    // loadData(paths, vHistoryRemainder, SolutionDataType::RemainderHistory);

    // Step 4: solve computational problem in Fourier domain using pseudo-spectral time-stepping method

    /* 
    classes:
    Solver	       methods: initialize/load optIC, load data, check CFL, Nonlinear dealisaing, IMEX RK4, save {fwdALL (optimizer)}, save {bwd,optTC (every time)}, save {energy,spectrum (last time)}, load {fwd (bwd)}, del {fwd (last time)}
    Optimizer	   methods: RCG, brent, save {optIC, branch, diags, brent}
    */ 

    // Solver solver(params);

    // solver.setInitialCondition(vStateInitial);
            // fill initial condition
        // vStateInitial.setInitialData(params, paths, vStateInitial)// init using str_guess, b_NumericalContinuation==0,d_OptimalTimeWindow
            // look for dir optIC or use InitialGuessName

    // solver.solvePDE(vStateInitial,vHistoryIntermediate,vHistoryRemainder);
        // solve forward PDE
        // computeSolution(params, paths, vStateInitial, vHistoryIntermediate, vHistoryRemainder) // then if b_opt == 0 updateDirData

    // solve optimization problem
    // if (params.b_OptimizeSolution == 1) {
    // Optimizer optimizer(params);
    // optimizer.maximizeEnergy(params, paths, solver,vStateInitial,vHistoryIntermediate,vHistoryRemainder);
        // optimizeSolution(params, paths, vStateInitial, vHistoryIntermediate, vHistoryRemainder) // then updateDirectoryData(params,const &paths)// save diags and delete forward solution in dir
    // }



    // Step 5: Clean up  
    // deleteData(paths);
    timer.stop();

    return 0;
}