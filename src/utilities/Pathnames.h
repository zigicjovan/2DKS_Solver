#ifndef PATHNAMES_H
#define PATHNAMES_H

#include "Parameters.h"
#include "MPIContext.h"

#include <filesystem>
#include <sstream>

using namespace std;

class Pathnames {
private:

    bool _useTempDir;
    bool _tempDirChoice;

    // directories 
    filesystem::path _dirData;
    filesystem::path _dirForwardSolution;
    filesystem::path _dirBackwardSolution;
    filesystem::path _dirFourierSpectrumEvolution;
    filesystem::path _dirEnergyEvolution;
    filesystem::path _dirInitialData;
    filesystem::path _dirTerminalData;
    filesystem::path _dirSolutionBranches;
    filesystem::path _dirOptimizationDiagnostics;
    filesystem::path _dirOptimizationLineSearch;

    filesystem::path _dirDataTemp;
    filesystem::path _dirForwardSolutionTemp;
    filesystem::path _dirBackwardSolutionTemp;
    filesystem::path _dirFourierSpectrumEvolutionTemp;
    filesystem::path _dirEnergyEvolutionTemp;
    filesystem::path _dirInitialDataTemp;
    filesystem::path _dirTerminalDataTemp;
    filesystem::path _dirSolutionBranchesTemp;
    filesystem::path _dirOptimizationDiagnosticsTemp;
    filesystem::path _dirOptimizationLineSearchTemp;

    // filenames
    ostringstream _strTestcase;
    ostringstream _strTestcaseGenericTime;
    ostringstream _strTestcaseBranch;
    ostringstream _strTestcaseInitialEnergyPowerLaw;
    ostringstream _strTestcaseDomainSizePowerLaw;
    ostringstream _strTestcaseEnergyTimeWindowPowerLaw;
    ostringstream _strTestcaseDomainTimeWindowPowerLaw;

    filesystem::path _fForwardSolution; // note: multiple files 
    filesystem::path _fBackwardSolution;
    filesystem::path _fFourierSpectrumEvolution;
    filesystem::path _fEnergyEvolution;
    filesystem::path _fInitialData;
    filesystem::path _fTerminalData;
    filesystem::path _fSolutionBranches; 
    filesystem::path _fInitialEnergyPowerLaw; 
    filesystem::path _fDomainSizePowerLaw; 
    filesystem::path _fEnergyTimeWindowPowerLaw;
    filesystem::path _fDomainTimeWindowPowerLaw; 
    filesystem::path _fOptimizationDiagnostics;
    filesystem::path _fOptimizationLineSearch;

    filesystem::path _fForwardSolutionTemp; // note: multiple files 
    filesystem::path _fBackwardSolutionTemp;
    filesystem::path _fFourierSpectrumEvolutionTemp;
    filesystem::path _fEnergyEvolutionTemp;
    filesystem::path _fInitialDataTemp;
    filesystem::path _fTerminalDataTemp;
    filesystem::path _fSolutionBranchesTemp; 
    filesystem::path _fInitialEnergyPowerLawTemp; 
    filesystem::path _fDomainSizePowerLawTemp; 
    filesystem::path _fEnergyTimeWindowPowerLawTemp;
    filesystem::path _fDomainTimeWindowPowerLawTemp; 
    filesystem::path _fOptimizationDiagnosticsTemp;
    filesystem::path _fOptimizationLineSearchTemp;


public:
    Pathnames(const Parameters &params, const MPIContext& mpi);

    void setTempDir(bool useTempDir);
    void resetTempDir();

    const filesystem::path& getDirData() const;
    const filesystem::path& getDirForwardSolution() const;
    const filesystem::path& getDirBackwardSolution() const;
    const filesystem::path& getDirFourierSpectrumEvolution() const;
    const filesystem::path& getDirEnergyEvolution() const;
    const filesystem::path& getDirInitialData() const;
    const filesystem::path& getDirTerminalData() const;
    const filesystem::path& getDirSolutionBranches() const;
    const filesystem::path& getDirOptimizationDiagnostics() const;
    const filesystem::path& getDirOptimizationLineSearch() const;

    string getTestcase() const;
    string getTestcaseGenericTime() const;
    string getTestcaseBranch() const;
    string getTestcaseInitialEnergyPowerLaw() const;
    string getTestcaseDomainSizePowerLaw() const;
    string getTestcaseEnergyTimeWindowPowerLaw() const;
    string getTestcaseDomainTimeWindowPowerLaw() const;

    const filesystem::path& getForwardSolutionFile() const;
    const filesystem::path& getBackwardSolutionFile() const;
    const filesystem::path& getFourierSpectrumEvolutionFile() const;
    const filesystem::path& getEnergyEvolutionFile() const;
    const filesystem::path& getInitialDataFile() const;
    void setInitialDataFile(const filesystem::path& initialDataFile);
    const filesystem::path& getTerminalDataFile() const;
    const filesystem::path& getSolutionBranchesFile() const;
    const filesystem::path& getInitialEnergyPowerLawFile() const;
    const filesystem::path& getDomainSizePowerLawFile() const;
    const filesystem::path& getEnergyTimeWindowPowerLawFile() const;
    const filesystem::path& getDomainTimeWindowPowerLawFile() const;
    const filesystem::path& getOptimizationDiagnosticsFile() const;
    const filesystem::path& getOptimizationLineSearchFile() const;
};

#endif  