// what: declare prototype for directory file and path names based on parameters
// what: data processed in C++, media (figures/movies) post-processed in MATLAB

#ifndef PATHNAMES_H  
#define PATHNAMES_H   

#include <filesystem>

struct Pathnames // how: only one file or directory per string unless noted below
{
    // directories 
    std::filesystem::path dir_Data;
    std::filesystem::path dir_ForwardSolution;
    std::filesystem::path dir_BackwardSolution;
    std::filesystem::path dir_FourierSpectrumEvolution;
    std::filesystem::path dir_EnergyEvolution;
    std::filesystem::path dir_OptimalInitialData;
    std::filesystem::path dir_OptimalTerminalData;
    std::filesystem::path dir_OptimalSolutionBranches;
    std::filesystem::path dir_OptimizationDiagnostics;
    std::filesystem::path dir_OptimizationLineSearch;

    // filenames
    std::filesystem::path file_ForwardSolution; // note: generates a sequence of files depending on problem size 
    std::filesystem::path file_BackwardSolution;
    std::filesystem::path file_FourierSpectrumEvolution;
    std::filesystem::path file_EnergyEvolution;
    std::filesystem::path file_OptimalInitialData;
    std::filesystem::path file_OptimalTerminalData;
    std::filesystem::path file_OptimalSolutionBranches;
    std::filesystem::path file_OptimizationDiagnostics;
    std::filesystem::path file_OptimizationLineSearch;
};

#endif  