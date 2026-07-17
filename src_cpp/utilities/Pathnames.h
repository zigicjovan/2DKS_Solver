// data processed in C++, media (figures/movies) post-processed in MATLAB

#ifndef PATHNAMES_H
#define PATHNAMES_H

#include "Parameters.h"

#include <filesystem>

using namespace std;

class Pathnames {
private:
public:
    Pathnames(const Parameters &params);

    // directories 
    filesystem::path dirData;
    filesystem::path dirForwardSolution;
    filesystem::path dirBackwardSolution;
    filesystem::path dirFourierSpectrumEvolution;
    filesystem::path dirEnergyEvolution;
    filesystem::path dirOptimalInitialData;
    filesystem::path dirOptimalTerminalData;
    filesystem::path dirOptimalSolutionBranch;
    filesystem::path dirOptimizationDiagnostics;
    filesystem::path dirOptimizationLineSearch;

    // filenames
    ostringstream strTestcase;
    ostringstream strTestcaseGenericTime;
    ostringstream strTestcaseBranch;
    ostringstream strTestcaseInitialEnergyPowerLaw;
    ostringstream strTestcaseDomainSizePowerLaw;
    ostringstream strTestcaseEnergyTimeWindowPowerLaw;
    ostringstream strTestcaseDomainTimeWindowPowerLaw;
    filesystem::path fForwardSolution; // note: multiple files 
    filesystem::path fBackwardSolution;
    filesystem::path fFourierSpectrumEvolution;
    filesystem::path fEnergyEvolution;
    filesystem::path fOptimalInitialData;
    filesystem::path fOptimalTerminalData;
    filesystem::path fOptimalSolutionBranch; 
    filesystem::path fOptimalSolutionInitialEnergyPowerLaw; 
    filesystem::path fOptimalSolutionDomainSizePowerLaw; 
    filesystem::path fOptimalSolutionEnergyTimeWindowPowerLaw;
    filesystem::path fOptimalSolutionDomainTimeWindowPowerLaw; 
    filesystem::path fOptimizationDiagnostics;
    filesystem::path fOptimizationLineSearch;
};

#endif  