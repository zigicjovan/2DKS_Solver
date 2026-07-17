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
    
    strTestcaseGenericTime << "_IC_" << params.strInitialGuessName 
        << "_N1_" << params.iGridSize1 
        << "_N2_" << params.iGridSize2     
        << "_dt_" << scientific << setprecision(1) << params.dTimeStep         
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2;

    strTestcaseBranch << "_IC_" << params.strInitialGuessName       
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy    
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2;

    strTestcaseInitialEnergyPowerLaw << "_IC_" << params.strInitialGuessName         
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2;

    strTestcaseDomainSizePowerLaw << "_IC_" << params.strInitialGuessName       
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy;

    strTestcaseEnergyTimeWindowPowerLaw << "_IC_" << params.strInitialGuessName       
        << "_ell1_" << fixed << setprecision(2) << params.dDomainFactor1     
        << "_ell2_" << fixed << setprecision(2) << params.dDomainFactor2;

    strTestcaseDomainTimeWindowPowerLaw << "_IC_" << params.strInitialGuessName       
        << "_K_" << scientific << setprecision(1) << params.dInitialEnergy;

    // create directories if they do not exist
    dirData = "Data";
    dirForwardSolution = dirData / "ForwardSolution" / strTestcase.str() ;
    dirBackwardSolution = dirData / "BackwardSolution";
    dirFourierSpectrumEvolution = dirData / "FourierSpectrumEvolution";
    dirEnergyEvolution = dirData / "EnergyEvolution";
    dirOptimalInitialData = dirData / "OptimalInitialData";
    dirOptimalTerminalData = dirData / "OptimalTerminalData";
    dirOptimalSolutionBranch = dirData / "OptimalSolutionBranch";
    dirOptimizationDiagnostics = dirData / "OptimizationDiagnostics";
    dirOptimizationLineSearch = dirData / "OptimizationLineSearch";

    filesystem::create_directories(dirData);
    filesystem::create_directories(dirForwardSolution);
    filesystem::create_directories(dirBackwardSolution);
    filesystem::create_directories(dirFourierSpectrumEvolution);
    filesystem::create_directories(dirEnergyEvolution);
    filesystem::create_directories(dirOptimalInitialData);
    filesystem::create_directories(dirOptimalTerminalData);
    filesystem::create_directories(dirOptimalSolutionBranch);
    filesystem::create_directories(dirOptimizationDiagnostics);
    filesystem::create_directories(dirOptimizationLineSearch);

    fForwardSolution = dirForwardSolution / ( "fwd" + strTestcase.str() ); 
    fBackwardSolution = dirBackwardSolution / ( "gradJ" + strTestcase.str() );
    fFourierSpectrumEvolution = dirFourierSpectrumEvolution / ( "spectrum" + strTestcase.str() );
    fEnergyEvolution = dirEnergyEvolution / ( "energy" + strTestcase.str() );
    fOptimalInitialData = dirOptimalInitialData / ( "fwdIC" + strTestcase.str() );
    fOptimalTerminalData = dirOptimalTerminalData / ( "fwdTC" + strTestcase.str() );
    fOptimalSolutionBranch = dirOptimalSolutionBranch / ( "branch" + strTestcaseBranch.str() );
    fOptimalSolutionInitialEnergyPowerLaw = dirOptimalSolutionBranch / ( "powerlawK" + strTestcaseInitialEnergyPowerLaw.str() );
    fOptimalSolutionDomainSizePowerLaw = dirOptimalSolutionBranch / ( "powerlawL" + strTestcaseDomainSizePowerLaw.str() );
    fOptimalSolutionEnergyTimeWindowPowerLaw = dirOptimalSolutionBranch / ( "powerlawTK" + strTestcaseEnergyTimeWindowPowerLaw.str() );
    fOptimalSolutionDomainTimeWindowPowerLaw = dirOptimalSolutionBranch / ( "powerlawTL" + strTestcaseDomainTimeWindowPowerLaw.str() );
    fOptimizationDiagnostics = dirOptimizationDiagnostics / ( "diagnostics" + strTestcase.str() );
    fOptimizationLineSearch = dirOptimizationLineSearch / ( "linesearch" + strTestcase.str() );

    cout << "Directory name: " << strTestcase.str() << endl;
}