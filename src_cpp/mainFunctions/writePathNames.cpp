#include "Pathnames.h"
#include "Parameters.h"
#include "writePathNames.h"

#include <iostream>  
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>

Pathnames writePathNames(const Parameters &params) {
    Pathnames paths;

    // prepare filenames
    paths.strTestcase << "_IC_" << params.strInitialGuessName 
        << "_N1_" << params.iGridSize1 
        << "_N2_" << params.iGridSize2     
        << "_dt_" << std::scientific << std::setprecision(1) << params.dTimeStep         
        << "_K_" << std::scientific << std::setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << std::fixed << std::setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << std::fixed << std::setprecision(2) << params.dDomainFactor2     
        << "_T_" << std::scientific << std::setprecision(2) << params.dTimeWindow       
        << "_opt_" << params.bOptimizeSolution       
        << "_tol_" << std::scientific << std::setprecision(0) << params.dOptimizationTolerance  
        << "_cont_" << params.bNumericalContinuation  
        << "_optT_" << std::scientific << std::setprecision(2) << params.dOptimalTimeWindow      
        << ".dat";
    
    paths.strTestcaseGeneric << "_IC_" << params.strInitialGuessName 
        << "_N1_" << params.iGridSize1 
        << "_N2_" << params.iGridSize2     
        << "_dt_" << std::scientific << std::setprecision(1) << params.dTimeStep         
        << "_K_" << std::scientific << std::setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << std::fixed << std::setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << std::fixed << std::setprecision(2) << params.dDomainFactor2;

    // create directories if they do not exist
    paths.dirData = "Data";
    paths.dirForwardSolution = paths.dirData / "ForwardSolution" / paths.strTestcase.str() ;
    paths.dirBackwardSolution = paths.dirData / "BackwardSolution";
    paths.dirFourierSpectrumEvolution = paths.dirData / "FourierSpectrumEvolution";
    paths.dirEnergyEvolution = paths.dirData / "EnergyEvolution";
    paths.dirOptimalInitialData = paths.dirData / "OptimalInitialData";
    paths.dirOptimalTerminalData = paths.dirData / "OptimalTerminalData";
    paths.dirOptimalSolutionBranches = paths.dirData / "OptimalSolutionBranches";
    paths.dirOptimizationDiagnostics = paths.dirData / "OptimizationDiagnostics";
    paths.dirOptimizationLineSearch = paths.dirData / "OptimizationLineSearch";

    std::filesystem::create_directories(paths.dirData);
    std::filesystem::create_directories(paths.dirForwardSolution);
    std::filesystem::create_directories(paths.dirBackwardSolution);
    std::filesystem::create_directories(paths.dirFourierSpectrumEvolution);
    std::filesystem::create_directories(paths.dirEnergyEvolution);
    std::filesystem::create_directories(paths.dirOptimalInitialData);
    std::filesystem::create_directories(paths.dirOptimalTerminalData);
    std::filesystem::create_directories(paths.dirOptimalSolutionBranches);
    std::filesystem::create_directories(paths.dirOptimizationDiagnostics);
    std::filesystem::create_directories(paths.dirOptimizationLineSearch);

    paths.fForwardSolution = paths.dirForwardSolution / ( "fwd" + paths.strTestcase.str() ); 
    paths.fBackwardSolution = paths.dirBackwardSolution / ( "bwd" + paths.strTestcase.str() );
    paths.fFourierSpectrumEvolution = paths.dirFourierSpectrumEvolution / ( "four" + paths.strTestcase.str() );
    paths.fEnergyEvolution = paths.dirEnergyEvolution / ( "energy" + paths.strTestcase.str() );
    paths.fOptimalInitialData = paths.dirOptimalInitialData / ( "optIC" + paths.strTestcase.str() );
    paths.fOptimalTerminalData = paths.dirOptimalTerminalData / ( "optTC" + paths.strTestcase.str() );
    paths.fOptimalSolutionBranches = paths.dirOptimalSolutionBranches / ( "branch" + paths.strTestcase.str() );
    paths.fOptimizationDiagnostics = paths.dirOptimizationDiagnostics / ( "diag" + paths.strTestcase.str() );
    paths.fOptimizationLineSearch = paths.dirOptimizationLineSearch / ( "brent" + paths.strTestcase.str() );

    std::cout << "Directory name: " << paths.strTestcase.str() << std::endl;
    return paths;
}