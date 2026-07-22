#ifndef PATHNAMES_H
#define PATHNAMES_H

#include "Parameters.h"
#include "MPIContext.h"

#include <filesystem>
#include <sstream>

using namespace std;

class Pathnames {
private:
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

public:
    Pathnames(const Parameters &params, const MPIContext& mpi);

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