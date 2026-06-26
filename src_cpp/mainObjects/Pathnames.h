// data processed in C++, media (figures/movies) post-processed in MATLAB

#ifndef PATHNAMES_H
#define PATHNAMES_H

#include <filesystem>

struct Pathnames {
    // directories 
    std::filesystem::path dirData;
    std::filesystem::path dirForwardSolution;
    std::filesystem::path dirBackwardSolution;
    std::filesystem::path dirFourierSpectrumEvolution;
    std::filesystem::path dirEnergyEvolution;
    std::filesystem::path dirOptimalInitialData;
    std::filesystem::path dirOptimalTerminalData;
    std::filesystem::path dirOptimalSolutionBranches;
    std::filesystem::path dirOptimizationDiagnostics;
    std::filesystem::path dirOptimizationLineSearch;

    // filenames
    std::filesystem::path fForwardSolution; // note: multiple files 
    std::filesystem::path fBackwardSolution;
    std::filesystem::path fFourierSpectrumEvolution;
    std::filesystem::path fEnergyEvolution;
    std::filesystem::path fOptimalInitialData;
    std::filesystem::path fOptimalTerminalData;
    std::filesystem::path fOptimalSolutionBranches;
    std::filesystem::path fOptimizationDiagnostics;
    std::filesystem::path fOptimizationLineSearch;
};

#endif  