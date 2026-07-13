#include "Pathnames.h"
#include "Parameters.h"

#include <iostream>  
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>

using namespace std;

Pathnames::Pathnames(const Parameters &params) {

    // prepare filenames
    strTestcase << "_IC_" << params.strInitialGuessName 
        << "_N1_" << params.iGridSize1 
        << "_N2_" << params.iGridSize2     
        << "_dt_" << scientific << setprecision(1) << params.dTimeStep         
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2     
        << "_T_" << scientific << setprecision(2) << params.dTimeWindow       
        << "_opt_" << params.bOptimizeSolution       
        << "_tol_" << scientific << setprecision(0) << params.dOptimizationTolerance  
        << "_cont_" << params.bNumericalContinuation  
        << "_optT_" << scientific << setprecision(2) << params.dOptimalTimeWindow      
        << ".dat";
    
    strTestcaseGeneric << "_IC_" << params.strInitialGuessName 
        << "_N1_" << params.iGridSize1 
        << "_N2_" << params.iGridSize2     
        << "_dt_" << scientific << setprecision(1) << params.dTimeStep         
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2;

    // create directories if they do not exist
    dirData = "Data";
    dirForwardSolution = dirData / "ForwardSolution" / strTestcase.str() ;
    dirBackwardSolution = dirData / "BackwardSolution";
    dirFourierSpectrumEvolution = dirData / "FourierSpectrumEvolution";
    dirEnergyEvolution = dirData / "EnergyEvolution";
    dirOptimalInitialData = dirData / "OptimalInitialData";
    dirOptimalTerminalData = dirData / "OptimalTerminalData";
    dirOptimalSolutionBranches = dirData / "OptimalSolutionBranches";
    dirOptimizationDiagnostics = dirData / "OptimizationDiagnostics";
    dirOptimizationLineSearch = dirData / "OptimizationLineSearch";

    filesystem::create_directories(dirData);
    filesystem::create_directories(dirForwardSolution);
    filesystem::create_directories(dirBackwardSolution);
    filesystem::create_directories(dirFourierSpectrumEvolution);
    filesystem::create_directories(dirEnergyEvolution);
    filesystem::create_directories(dirOptimalInitialData);
    filesystem::create_directories(dirOptimalTerminalData);
    filesystem::create_directories(dirOptimalSolutionBranches);
    filesystem::create_directories(dirOptimizationDiagnostics);
    filesystem::create_directories(dirOptimizationLineSearch);

    fForwardSolution = dirForwardSolution / ( "fwd" + strTestcase.str() ); 
    fBackwardSolution = dirBackwardSolution / ( "bwd" + strTestcase.str() );
    fFourierSpectrumEvolution = dirFourierSpectrumEvolution / ( "four" + strTestcase.str() );
    fEnergyEvolution = dirEnergyEvolution / ( "energy" + strTestcase.str() );
    fOptimalInitialData = dirOptimalInitialData / ( "optIC" + strTestcase.str() );
    fOptimalTerminalData = dirOptimalTerminalData / ( "optTC" + strTestcase.str() );
    fOptimalSolutionBranches = dirOptimalSolutionBranches / ( "branch" + strTestcase.str() );
    fOptimizationDiagnostics = dirOptimizationDiagnostics / ( "diag" + strTestcase.str() );
    fOptimizationLineSearch = dirOptimizationLineSearch / ( "brent" + strTestcase.str() );

    cout << "Directory name: " << strTestcase.str() << endl;
}