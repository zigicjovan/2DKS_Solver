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
    filesystem::path dirOptimalSolutionBranches;
    filesystem::path dirOptimizationDiagnostics;
    filesystem::path dirOptimizationLineSearch;

    // filenames
    ostringstream strTestcase;
    ostringstream strTestcaseGeneric;
    filesystem::path fForwardSolution; // note: multiple files 
    filesystem::path fBackwardSolution;
    filesystem::path fFourierSpectrumEvolution;
    filesystem::path fEnergyEvolution;
    filesystem::path fOptimalInitialData;
    filesystem::path fOptimalTerminalData;
    filesystem::path fOptimalSolutionBranches;
    filesystem::path fOptimizationDiagnostics;
    filesystem::path fOptimizationLineSearch;
};

#endif  