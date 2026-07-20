#include "Parameters.h"
#include "Pathnames.h"
#include "FFTWPlanner.h"
#include "Timer.h"
#include "Solver.h"
#include "SolutionData.h"
#include <mpi.h>
#include <fftw3-mpi.h>
#include "MPIContext.h"

#include <iostream>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {

    // Step 1: allocate distributed computing resources
    MPI_Init(&argc, &argv);
    fftw_mpi_init();
    int mpiRank;
    int mpiSize;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
    
    {
        // Step 2: prepare test case
        Parameters params(argc, argv);
        MPIContext mpi(params.getGridSize1(), params.getGridSize2());
        params.getMathematicalOperators(mpi);
        Pathnames paths(params, mpi);
        FFTWPlanner fftwPlan(params, mpi);
        mpi.printDecomposition();

        // Step 3: track computation time
        Timer timer(mpi);
        timer.start(); 

        // Step 4: allocate storage for states and time-dependent solutions
        SolutionData vStateInitial(params, paths, mpi, InitialState);      
        vStateInitial.validateMPIConfiguration();
        SolutionData vStateTerminal(params, paths, mpi, TerminalState); 
        SolutionData vObjectiveGradient(params, paths, mpi, BackwardInitialState); 
        SolutionData vHistoryIntermediate(params, paths, mpi, IntermediateHistory);
        SolutionData vHistoryRemainder(params, paths, mpi, RemainderHistory);

        // Step 5: solve computational problem in Fourier domain using pseudo-spectral time-stepping method
        Solver solver(params, paths, fftwPlan, timer, mpi);
        solver.setSolutionState(SolveInitialState, vStateInitial);
        solver.setSolutionInTime(SolveForwardInTime, vStateInitial, vHistoryIntermediate, vHistoryRemainder, vStateTerminal);
        timer.printInterval("Forward problem initialized at ");
        solver.getOptimalSolution(OptimizeEnergyAmplification, vObjectiveGradient, vStateInitial, vHistoryIntermediate, vHistoryRemainder, vStateTerminal);  

        // Step 6: Stop timer and complete run  
        timer.stop();
    }

    // Step 7: deallocate distributed computing resources 
    fftw_mpi_cleanup();
    MPI_Finalize();

    return 0;
}