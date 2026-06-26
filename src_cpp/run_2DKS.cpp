// what: add paths to headers of supporting files
#include "readParameterSettings.h"
#include "writePathNames.h"
#include "Timer.h"
#include "SolutionData.h"
#include "updateDirectory.h"
// #include "computeSolution.h"
// #include "optimizeSolution.h"

int main(int argc, char* argv[])  
{
    // Step 1: prepare test case
    Parameters params = readParameterSettings(argc, argv);
    Pathnames paths = writePathNames(params);

    // Step 2: track computation time
    Timer timer;
    timer.start(); 
    // use timer.printInterval("");

    // Step 3: allocate storage for initial state and time-dependent solution
    SolutionData state_Initial(params, SolutionDataType::InitialState); 
    SolutionData history_Intermediate(params, SolutionDataType::IntermediateHistory); // why: if ( params.n_numericalSteps() > params.n_numericalStepsPerFile() )
    SolutionData history_Remainder(params, SolutionDataType::RemainderHistory); 
    /* use
    saveData(paths, state_Initial, SolutionDataType::InitialState);
    saveData(paths, history_Intermediate, SolutionDataType::IntermediateHistory, 1000);
    saveData(paths, history_Remainder, SolutionDataType::RemainderHistory, params.d_TimeWindow);
    loadData(paths, state_Initial, SolutionDataType::InitialState);
    loadData(paths, history_Intermediate, SolutionDataType::IntermediateHistory);
    loadData(paths, history_Remainder, SolutionDataType::RemainderHistory);
    */

    // Step 4: solve computational problem in Fourier domain using pseudo-spectral time-stepping method

    // fill initial condition
        // state_Initial.setInitialData(params,const &paths,&state_Initial)// init using str_guess, b_NumericalContinuation==0,d_OptimalTimeWindow
            // look for dir optIC or use InitialGuessName
    // solve forward PDE
        // computeSolution(params,const &paths,&state_Initial,&history_Intermediate,&history_Remaining) // then if b_opt == 0 updateDirData
    // solve optimization problem
        // if (params.b_OptimizeSolution == 1), 
        // optimizeSolution(params,const &paths,&state_Initial,&history_Intermediate,&history_Remaining) // then updateDirectoryData(params,const &paths)// save diags and delete forward solution in dir
    
    /* 
    classes:
    SolutionData   methods:	initialize/load optIC, load data 
    Solver	       methods: check CFL, Nonlinear dealisaing, IMEX RK4, save fwdALL (optimizer), save bwd,optTC (every time), save energy/spectrum (last time), load fwd (bwd), del fwd (last time)
    Optimizer	   methods: RCG, brent, save optIC,branch,diags,brent
    */ 

    // Step 5: Clean up  
    // deleteData(paths);
    timer.stop(); // print computation summary

    return 0;
}