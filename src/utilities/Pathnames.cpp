#include "Pathnames.h"
#include "Parameters.h"

#include <iostream>  
#include <string>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <cstdlib>

using namespace std;

Pathnames::Pathnames(const Parameters &params, const MPIContext& mpi) {

    // prepare filenames
    const string fType = ".dat";
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
        << "_rank_" << mpi.getSize();
    
    _strTestcaseGenericTime << "_IC_" << params.getInitialGuessName() 
        << "_N1_" << params.getGridSize1() 
        << "_N2_" << params.getGridSize2()     
        << "_dt_" << scientific << setprecision(1) << params.getTimeStep()         
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy()    
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2();

    _strTestcaseBranch << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy()    
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2();

    _strTestcaseInitialEnergyPowerLaw << "_IC_" << params.getInitialGuessName()          
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2();

    _strTestcaseDomainSizePowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy();

    _strTestcaseEnergyTimeWindowPowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_ell1_" << fixed << setprecision(2) << params.getDomainFactor1()      
        << "_ell2_" << fixed << setprecision(2) << params.getDomainFactor2();

    _strTestcaseDomainTimeWindowPowerLaw << "_IC_" << params.getInitialGuessName()        
        << "_K_" << scientific << setprecision(1) << params.getInitialEnergy();

    const char* scratch = getenv("SCRATCH");
    const filesystem::path dataRoot = scratch ? filesystem::path(scratch) / "Data" : filesystem::path("Data");
    const char* slurmTemp = getenv("SLURM_TMPDIR");
    const bool haveTmpDir = (slurmTemp != nullptr);
    filesystem::path dataRootTemp = dataRoot;
    _useTempDir = false;
    _tempDirChoice = _useTempDir;

    if (params.getOptimizeSolution() && haveTmpDir) {
        dataRootTemp = filesystem::path(slurmTemp) / "Data";

        const auto spaceInfo = filesystem::space(slurmTemp);
        const double availableGB = static_cast<double>(spaceInfo.available) / (1024.0 * 1024.0 * 1024.0);

        _useTempDir = (params.getCheckpointMemory() <= 0.8 * availableGB);
        _tempDirChoice = _useTempDir;
    }

    const string testcase = _strTestcase.str();

    // create directories if they do not exist
    _dirData = dataRoot / testcase ;
    _dirForwardSolution = _dirData / "ForwardSolution" ;
    _dirBackwardSolution = _dirData / "BackwardSolution";
    _dirFourierSpectrumEvolution = _dirData / "FourierSpectrumEvolution";
    _dirEnergyEvolution = _dirData / "EnergyEvolution";
    _dirInitialData = _dirData / "InitialData";
    _dirTerminalData = _dirData / "TerminalData";
    _dirSolutionBranches = dataRoot / "SolutionBranches" / testcase ;
    _dirOptimizationDiagnostics = _dirData / "OptimizationDiagnostics";
    _dirOptimizationLineSearch = _dirData / "OptimizationLineSearch";

    _fForwardSolution = _dirForwardSolution / ( "fwd" + fType ); 
    _fBackwardSolution = _dirBackwardSolution / ( "gradJ" + fType );
    _fFourierSpectrumEvolution = _dirFourierSpectrumEvolution / ( "spectrum" + fType );
    _fEnergyEvolution = _dirEnergyEvolution / ( "energy" + fType );
    _fInitialData = _dirInitialData / ( "fwdIC" + fType );
    _fTerminalData = _dirTerminalData / ( "fwdTC" + fType );
    _fSolutionBranches = _dirSolutionBranches / ( "branch" + _strTestcaseBranch.str() + fType );
    _fInitialEnergyPowerLaw = _dirSolutionBranches / ( "powerlawK" + _strTestcaseInitialEnergyPowerLaw.str() + fType );
    _fDomainSizePowerLaw = _dirSolutionBranches / ( "powerlawL" + _strTestcaseDomainSizePowerLaw.str() + fType );
    _fEnergyTimeWindowPowerLaw = _dirSolutionBranches / ( "powerlawTK" + _strTestcaseEnergyTimeWindowPowerLaw.str() + fType );
    _fDomainTimeWindowPowerLaw = _dirSolutionBranches / ( "powerlawTL" + _strTestcaseDomainTimeWindowPowerLaw.str() + fType );
    _fOptimizationDiagnostics = _dirOptimizationDiagnostics / ( "diagnostics" + fType );
    _fOptimizationLineSearch = _dirOptimizationLineSearch / ( "linesearch" + fType );

    // create directories if they do not exist
    _dirDataTemp = dataRootTemp / testcase ;
    _dirForwardSolutionTemp = _dirDataTemp / "ForwardSolution" ;
    _dirBackwardSolutionTemp = _dirDataTemp / "BackwardSolution";
    _dirFourierSpectrumEvolutionTemp = _dirDataTemp / "FourierSpectrumEvolution";
    _dirEnergyEvolutionTemp = _dirDataTemp / "EnergyEvolution";
    _dirInitialDataTemp = _dirDataTemp / "InitialData";
    _dirTerminalDataTemp = _dirDataTemp / "TerminalData";
    _dirSolutionBranchesTemp = dataRootTemp / "SolutionBranches" / testcase ;
    _dirOptimizationDiagnosticsTemp = _dirDataTemp / "OptimizationDiagnostics";
    _dirOptimizationLineSearchTemp = _dirDataTemp / "OptimizationLineSearch";

    _fForwardSolutionTemp = _dirForwardSolutionTemp / ( "fwd" + fType ); 
    _fBackwardSolutionTemp = _dirBackwardSolutionTemp / ( "gradJ" + fType );
    _fFourierSpectrumEvolutionTemp = _dirFourierSpectrumEvolutionTemp / ( "spectrum" + fType );
    _fEnergyEvolutionTemp = _dirEnergyEvolutionTemp / ( "energy" + fType );
    _fInitialDataTemp = _dirInitialDataTemp / ( "fwdIC" + fType );
    _fTerminalDataTemp = _dirTerminalDataTemp / ( "fwdTC" + fType );
    _fSolutionBranchesTemp = _dirSolutionBranchesTemp / ( "branch" + _strTestcaseBranch.str() + fType );
    _fInitialEnergyPowerLawTemp = _dirSolutionBranchesTemp / ( "powerlawK" + _strTestcaseInitialEnergyPowerLaw.str() + fType );
    _fDomainSizePowerLawTemp = _dirSolutionBranchesTemp / ( "powerlawL" + _strTestcaseDomainSizePowerLaw.str() + fType );
    _fEnergyTimeWindowPowerLawTemp = _dirSolutionBranchesTemp / ( "powerlawTK" + _strTestcaseEnergyTimeWindowPowerLaw.str() + fType );
    _fDomainTimeWindowPowerLawTemp = _dirSolutionBranchesTemp / ( "powerlawTL" + _strTestcaseDomainTimeWindowPowerLaw.str() + fType );
    _fOptimizationDiagnosticsTemp = _dirOptimizationDiagnosticsTemp / ( "diagnostics" + fType );
    _fOptimizationLineSearchTemp = _dirOptimizationLineSearchTemp / ( "linesearch" + fType );

    if (mpi.isRoot()) {
        filesystem::create_directories(_dirForwardSolution);
        filesystem::create_directories(_dirBackwardSolution);
        filesystem::create_directories(_dirFourierSpectrumEvolution);
        filesystem::create_directories(_dirEnergyEvolution);
        filesystem::create_directories(_dirInitialData);
        filesystem::create_directories(_dirTerminalData);
        filesystem::create_directories(_dirSolutionBranches);
        filesystem::create_directories(_dirOptimizationDiagnostics);
        filesystem::create_directories(_dirOptimizationLineSearch);

        filesystem::create_directories(_dirForwardSolutionTemp);
        filesystem::create_directories(_dirBackwardSolutionTemp);
        filesystem::create_directories(_dirFourierSpectrumEvolutionTemp);
        filesystem::create_directories(_dirEnergyEvolutionTemp);
        filesystem::create_directories(_dirInitialDataTemp);
        filesystem::create_directories(_dirTerminalDataTemp);
        filesystem::create_directories(_dirSolutionBranchesTemp);
        filesystem::create_directories(_dirOptimizationDiagnosticsTemp);
        filesystem::create_directories(_dirOptimizationLineSearchTemp);

        cout << "Writing temporary data to " << (_useTempDir ? "SLURM_TMPDIR" : "SCRATCH") << endl;
        cout << "Directory name: " << testcase << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void Pathnames::setTempDir(bool useTempDir) {
    _useTempDir = useTempDir;
}

void Pathnames::resetTempDir() {
    _useTempDir = _tempDirChoice;
}

const filesystem::path& Pathnames::getDirData() const {
    if (!_useTempDir) {
        return _dirData;
    }
    else {
        return _dirDataTemp;
    }
}

const filesystem::path& Pathnames::getDirForwardSolution() const {
    if (!_useTempDir) {
        return _dirForwardSolution;
    }
    else {
        return _dirForwardSolutionTemp;
    }
}

const filesystem::path& Pathnames::getDirBackwardSolution() const {
    if (!_useTempDir) {
        return _dirBackwardSolution;
    }
    else {
        return _dirBackwardSolutionTemp;
    }
}

const filesystem::path& Pathnames::getDirFourierSpectrumEvolution() const {
    if (!_useTempDir) {
        return _dirFourierSpectrumEvolution;
    }
    else {
        return _dirFourierSpectrumEvolutionTemp;
    }
}

const filesystem::path& Pathnames::getDirEnergyEvolution() const {
    if (!_useTempDir) {
        return _dirEnergyEvolution;
    }
    else {
        return _dirEnergyEvolutionTemp;
    }
}

const filesystem::path& Pathnames::getDirInitialData() const {
    if (!_useTempDir) {
        return _dirInitialData;
    }
    else {
        return _dirInitialDataTemp;
    }
}

const filesystem::path& Pathnames::getDirTerminalData() const {
    if (!_useTempDir) {
        return _dirTerminalData;
    }
    else {
        return _dirTerminalDataTemp;
    }
}

const filesystem::path& Pathnames::getDirSolutionBranches() const {
    if (!_useTempDir) {
        return _dirSolutionBranches;
    }
    else {
        return _dirSolutionBranchesTemp;
    }
}

const filesystem::path& Pathnames::getDirOptimizationDiagnostics() const {
    if (!_useTempDir) {
        return _dirOptimizationDiagnostics;
    }
    else {
        return _dirOptimizationDiagnosticsTemp;
    }
}

const filesystem::path& Pathnames::getDirOptimizationLineSearch() const {
    if (!_useTempDir) {
        return _dirOptimizationLineSearch;
    }
    else {
        return _dirOptimizationLineSearchTemp;
    }
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
    if (!_useTempDir) {
        return _fForwardSolution;
    }
    else {
        return _fForwardSolutionTemp;
    }
}

const filesystem::path& Pathnames::getBackwardSolutionFile() const {
    if (!_useTempDir) {
        return _fBackwardSolution;
    }
    else {
        return _fBackwardSolutionTemp;
    }
}

const filesystem::path& Pathnames::getFourierSpectrumEvolutionFile() const {
    if (!_useTempDir) {
        return _fFourierSpectrumEvolution;
    }
    else {
        return _fFourierSpectrumEvolutionTemp;
    }
}

const filesystem::path& Pathnames::getEnergyEvolutionFile() const {
    if (!_useTempDir) {
        return _fEnergyEvolution;
    }
    else {
        return _fEnergyEvolutionTemp;
    }
}

const filesystem::path& Pathnames::getInitialDataFile() const {
    if (!_useTempDir) {
        return _fInitialData;
    }
    else {
        return _fInitialDataTemp;
    }
}

void Pathnames::setInitialDataFile(const filesystem::path& initialDataFile) {
    _fInitialData = initialDataFile;
}

const filesystem::path& Pathnames::getTerminalDataFile() const {
    if (!_useTempDir) {
        return _fTerminalData;
    }
    else {
        return _fTerminalDataTemp;
    }
}

const filesystem::path& Pathnames::getSolutionBranchesFile() const {
    if (!_useTempDir) {
        return _fSolutionBranches;
    }
    else {
        return _fSolutionBranchesTemp;
    }
}

const filesystem::path& Pathnames::getInitialEnergyPowerLawFile() const {
    if (!_useTempDir) {
        return _fInitialEnergyPowerLaw;
    }
    else {
        return _fInitialEnergyPowerLawTemp;
    }
}

const filesystem::path& Pathnames::getDomainSizePowerLawFile() const {
    if (!_useTempDir) {
        return _fDomainSizePowerLaw;
    }
    else {
        return _fDomainSizePowerLawTemp;
    }
}

const filesystem::path& Pathnames::getEnergyTimeWindowPowerLawFile() const {
    if (!_useTempDir) {
        return _fEnergyTimeWindowPowerLaw;
    }
    else {
        return _fEnergyTimeWindowPowerLawTemp;
    }
}

const filesystem::path& Pathnames::getDomainTimeWindowPowerLawFile() const {
    if (!_useTempDir) {
        return _fDomainTimeWindowPowerLaw;
    }
    else {
        return _fDomainTimeWindowPowerLawTemp;
    }
}

const filesystem::path& Pathnames::getOptimizationDiagnosticsFile() const {
    if (!_useTempDir) {
        return _fOptimizationDiagnostics;
    }
    else {
        return _fOptimizationDiagnosticsTemp;
    }
}

const filesystem::path& Pathnames::getOptimizationLineSearchFile() const {
    if (!_useTempDir) {
        return _fOptimizationLineSearch;
    }
    else {
        return _fOptimizationLineSearchTemp;
    }
}