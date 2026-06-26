#include "Pathnames.h"
#include "Parameters.h"
#include "writePathNames.h"

#include <iostream>  
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>

Pathnames writePathNames(const Parameters &params)
{
    Pathnames paths;

    // prepare filenames
    std::ostringstream testcase;
    testcase    << "_IC_"      << params.str_InitalGuessName 
                << "_N1_"      << params.n_GridSize_1 
                << "_N2_"      << params.n_GridSize_2     
                << "_dt_"        << std::scientific << std::setprecision(1) << params.d_TimeStep         
                << "_K_"         << std::scientific << std::setprecision(1) << params.d_InitialEnergy    
                << "_ell1_"      << std::fixed      << std::setprecision(2) << params.d_DomainFactor_1     
                << "_ell2_"      << std::fixed      << std::setprecision(2) << params.d_DomainFactor_2     
                << "_T_"         << std::scientific << std::setprecision(2) << params.d_TimeWindow       
                << "_opt_"       << params.b_OptimizeSolution       
                << "_tol_"       << std::scientific << std::setprecision(0) << params.d_OptimizationTolerance  
                << "_cont_"      << params.b_NumericalContinuation  
                << "_optT_"      << std::scientific << std::setprecision(2) << params.d_OptimalTimeWindow      
                << ".dat";

    // create directories if they do not exist
    paths.dir_Data                      = "Data";
    paths.dir_ForwardSolution           = paths.dir_Data / "ForwardSolution" / testcase.str() ;
    paths.dir_BackwardSolution          = paths.dir_Data / "BackwardSolution";
    paths.dir_FourierSpectrumEvolution  = paths.dir_Data / "FourierSpectrumEvolution";
    paths.dir_EnergyEvolution           = paths.dir_Data / "EnergyEvolution";
    paths.dir_OptimalInitialData        = paths.dir_Data / "OptimalInitialData";
    paths.dir_OptimalTerminalData       = paths.dir_Data / "OptimalTerminalData";
    paths.dir_OptimalSolutionBranches   = paths.dir_Data / "OptimalSolutionBranches";
    paths.dir_OptimizationDiagnostics   = paths.dir_Data / "OptimizationDiagnostics";
    paths.dir_OptimizationLineSearch    = paths.dir_Data / "OptimizationLineSearch";


    std::filesystem::create_directories(paths.dir_Data);
    std::filesystem::create_directories(paths.dir_ForwardSolution);
    std::filesystem::create_directories(paths.dir_BackwardSolution);
    std::filesystem::create_directories(paths.dir_FourierSpectrumEvolution);
    std::filesystem::create_directories(paths.dir_EnergyEvolution);
    std::filesystem::create_directories(paths.dir_OptimalInitialData);
    std::filesystem::create_directories(paths.dir_OptimalTerminalData);
    std::filesystem::create_directories(paths.dir_OptimalSolutionBranches);
    std::filesystem::create_directories(paths.dir_OptimizationDiagnostics);
    std::filesystem::create_directories(paths.dir_OptimizationLineSearch);

    paths.file_ForwardSolution          = paths.dir_ForwardSolution / ( "fwd" + testcase.str() ); 
    paths.file_BackwardSolution         = paths.dir_BackwardSolution / ( "bwd" + testcase.str() );
    paths.file_FourierSpectrumEvolution = paths.dir_FourierSpectrumEvolution / ( "four" + testcase.str() );
    paths.file_EnergyEvolution          = paths.dir_EnergyEvolution / ( "energy" + testcase.str() );
    paths.file_OptimalInitialData       = paths.dir_OptimalInitialData / ( "optIC" + testcase.str() );
    paths.file_OptimalTerminalData      = paths.dir_OptimalTerminalData / ( "optTC" + testcase.str() );
    paths.file_OptimalSolutionBranches  = paths.dir_OptimalSolutionBranches / ( "branch" + testcase.str() );
    paths.file_OptimizationDiagnostics  = paths.dir_OptimizationDiagnostics / ( "diag" + testcase.str() );
    paths.file_OptimizationLineSearch   = paths.dir_OptimizationLineSearch / ( "brent" + testcase.str() );

    std::cout << "Testcase:\n" << testcase.str() << std::endl;

    return paths;
}