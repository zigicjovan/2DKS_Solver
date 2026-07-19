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
    _strTestcase << "_IC_" << params.getInitialGuessName() 
        << "_N1_" << params.getGridSize1() 
        << "_N2_" << params.getGridSize2()     
        << "_dt_" << scientific << setprecision(1) << params.getTimeStep()         
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy()    
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()     
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2()     
        << "_T_" << scientific << setprecision(2) << params.getTimeWindow()       
        << "_opt_" << params.getOptimizeSolution()       
        << "_tol_" << scientific << setprecision(0) << params.getOptimizationTolerance()  
        << "_cont_" << params.getNumericalContinuation()  
        << "_optT_" << scientific << setprecision(2) << params.getOptimalTimeWindow()      
        << ".dat";
    
    _strTestcaseGenericTime << "_IC_" << params.getInitialGuessName() 
        << "_N1_" << params.getGridSize1() 
        << "_N2_" << params.getGridSize2()     
        << "_dt_" << scientific << setprecision(1) << params.getTimeStep()         
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy()    
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2() ;

    _strTestcaseBranch << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy()    
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2()  << ".dat";

    _strTestcaseInitialEnergyPowerLaw << "_IC_" << params.getInitialGuessName()          
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2()  << ".dat";

    _strTestcaseDomainSizePowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy() << ".dat";

    _strTestcaseEnergyTimeWindowPowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2()  << ".dat";

    _strTestcaseDomainTimeWindowPowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy() << ".dat";

    // create directories if they do not exist
    _dirData = "Data";
    _dirForwardSolution = _dirData / "ForwardSolution" / _strTestcase.str() ;
    _dirBackwardSolution = _dirData / "BackwardSolution";
    _dirFourierSpectrumEvolution = _dirData / "FourierSpectrumEvolution";
    _dirEnergyEvolution = _dirData / "EnergyEvolution";
    _dirInitialData = _dirData / "InitialData";
    _dirTerminalData = _dirData / "TerminalData";
    _dirSolutionBranches = _dirData / "SolutionBranches" / _strTestcase.str() ;
    _dirOptimizationDiagnostics = _dirData / "OptimizationDiagnostics";
    _dirOptimizationLineSearch = _dirData / "OptimizationLineSearch";

    filesystem::create_directories(_dirData);
    filesystem::create_directories(_dirForwardSolution);
    filesystem::create_directories(_dirBackwardSolution);
    filesystem::create_directories(_dirFourierSpectrumEvolution);
    filesystem::create_directories(_dirEnergyEvolution);
    filesystem::create_directories(_dirInitialData);
    filesystem::create_directories(_dirTerminalData);
    filesystem::create_directories(_dirSolutionBranches);
    filesystem::create_directories(_dirOptimizationDiagnostics);
    filesystem::create_directories(_dirOptimizationLineSearch);

    _fForwardSolution = _dirForwardSolution / ( "fwd" + _strTestcase.str() ); 
    _fBackwardSolution = _dirBackwardSolution / ( "gradJ" + _strTestcase.str() );
    _fFourierSpectrumEvolution = _dirFourierSpectrumEvolution / ( "spectrum" + _strTestcase.str() );
    _fEnergyEvolution = _dirEnergyEvolution / ( "energy" + _strTestcase.str() );
    _fInitialData = _dirInitialData / ( "fwdIC" + _strTestcase.str() );
    _fTerminalData = _dirTerminalData / ( "fwdTC" + _strTestcase.str() );
    _fSolutionBranches = _dirSolutionBranches / ( "branch" + _strTestcaseBranch.str() );
    _fInitialEnergyPowerLaw = _dirSolutionBranches / ( "powerlawK" + _strTestcaseInitialEnergyPowerLaw.str() );
    _fDomainSizePowerLaw = _dirSolutionBranches / ( "powerlawL" + _strTestcaseDomainSizePowerLaw.str() );
    _fEnergyTimeWindowPowerLaw = _dirSolutionBranches / ( "powerlawTK" + _strTestcaseEnergyTimeWindowPowerLaw.str() );
    _fDomainTimeWindowPowerLaw = _dirSolutionBranches / ( "powerlawTL" + _strTestcaseDomainTimeWindowPowerLaw.str() );
    _fOptimizationDiagnostics = _dirOptimizationDiagnostics / ( "diagnostics" + _strTestcase.str() );
    _fOptimizationLineSearch = _dirOptimizationLineSearch / ( "linesearch" + _strTestcase.str() );

    cout << "Directory name: " << _strTestcase.str() << endl;
}

const filesystem::path& Pathnames::getDirData() const {
    return _dirData;
}

const filesystem::path& Pathnames::getDirForwardSolution() const {
    return _dirForwardSolution;
}

const filesystem::path& Pathnames::getDirBackwardSolution() const {
    return _dirBackwardSolution;
}

const filesystem::path& Pathnames::getDirFourierSpectrumEvolution() const {
    return _dirFourierSpectrumEvolution;
}

const filesystem::path& Pathnames::getDirEnergyEvolution() const {
    return _dirEnergyEvolution;
}

const filesystem::path& Pathnames::getDirInitialData() const {
    return _dirInitialData;
}

const filesystem::path& Pathnames::getDirTerminalData() const {
    return _dirTerminalData;
}

const filesystem::path& Pathnames::getDirSolutionBranches() const {
    return _dirSolutionBranches;
}

const filesystem::path& Pathnames::getDirOptimizationDiagnostics() const {
    return _dirOptimizationDiagnostics;
}

const filesystem::path& Pathnames::getDirOptimizationLineSearch() const {
    return _dirOptimizationLineSearch;
}

string Pathnames::getTestcase() const {
    return _strTestcase.str();
}

string Pathnames::getTestcaseGenericTime() const {
    return _strTestcaseGenericTime.str();
}

string Pathnames::getTestcaseBranch() const {
    return _strTestcaseBranch.str();
}

string Pathnames::getTestcaseInitialEnergyPowerLaw() const {
    return _strTestcaseInitialEnergyPowerLaw.str();
}

string Pathnames::getTestcaseDomainSizePowerLaw() const {
    return _strTestcaseDomainSizePowerLaw.str();
}

string Pathnames::getTestcaseEnergyTimeWindowPowerLaw() const {
    return _strTestcaseEnergyTimeWindowPowerLaw.str();
}

string Pathnames::getTestcaseDomainTimeWindowPowerLaw() const {
    return _strTestcaseDomainTimeWindowPowerLaw.str();
}

const filesystem::path& Pathnames::getForwardSolutionFile() const {
    return _fForwardSolution;
}

const filesystem::path& Pathnames::getBackwardSolutionFile() const {
    return _fBackwardSolution;
}

const filesystem::path& Pathnames::getFourierSpectrumEvolutionFile() const {
    return _fFourierSpectrumEvolution;
}

const filesystem::path& Pathnames::getEnergyEvolutionFile() const {
    return _fEnergyEvolution;
}

const filesystem::path& Pathnames::getInitialDataFile() const {
    return _fInitialData;
}

void Pathnames::setInitialDataFile(const filesystem::path& initialDataFile) {
    _fInitialData = initialDataFile;
}

const filesystem::path& Pathnames::getTerminalDataFile() const {
    return _fTerminalData;
}

const filesystem::path& Pathnames::getSolutionBranchesFile() const {
    return _fSolutionBranches;
}

const filesystem::path& Pathnames::getInitialEnergyPowerLawFile() const {
    return _fInitialEnergyPowerLaw;
}

const filesystem::path& Pathnames::getDomainSizePowerLawFile() const {
    return _fDomainSizePowerLaw;
}

const filesystem::path& Pathnames::getEnergyTimeWindowPowerLawFile() const {
    return _fEnergyTimeWindowPowerLaw;
}

const filesystem::path& Pathnames::getDomainTimeWindowPowerLawFile() const {
    return _fDomainTimeWindowPowerLaw;
}

const filesystem::path& Pathnames::getOptimizationDiagnosticsFile() const {
    return _fOptimizationDiagnostics;
}

const filesystem::path& Pathnames::getOptimizationLineSearchFile() const {
    return _fOptimizationLineSearch;
}