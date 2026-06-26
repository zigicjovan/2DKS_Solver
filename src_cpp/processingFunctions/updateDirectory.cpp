#include "Pathnames.h"
#include "SolutionData.h"

#include <filesystem>
#include <fstream>
#include <iostream> 

std::filesystem::path appendTimeStep(const std::filesystem::path& path, double dCurrentT) {
    return path.parent_path() / (path.stem().string() + "_" + std::to_string(dCurrentT) + path.extension().string());
}

void loadData(const Pathnames& paths, SolutionData& storedData, SolutionDataType storedDataType, double dCurrentT) {   
    /*
    SolutionDataTypes:
    // Global scope, store as binary:
    InitialState,               // why: load for fwd and IC
    TerminalState,              // why: load for bwd
    BackwardInitialState,       // why: load for optimization (i.e. objective functional gradient)
    IntermediateHistory,        // why: load for bwd
    RemainderHistory,           // why: load for bwd
    
    // Local scope, do this inline:
    SolutionBranch,             // why: load for optimization
    */  

    switch (storedDataType) {
        case InitialState: {
            std::ifstream file(paths.fOptimalInitialData, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex));// why: save as binary
            break;
        }
        case TerminalState: {
            std::ifstream file(paths.fOptimalTerminalData, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case BackwardInitialState: {
            std::ifstream file(paths.fBackwardSolution, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case IntermediateHistory:   
        case RemainderHistory: {
            std::ifstream file(appendTimeStep(paths.fForwardSolution, dCurrentT), std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
    }
}

void saveData(const Pathnames& paths, SolutionData& storedData, SolutionDataType storedDataType, double dCurrentT) {   
    /*
    SolutionDataTypes:
    // Global scope, store as binary:
    InitialState,               // why: save optIC
    TerminalState,              // why: save optTC
    BackwardInitialState,       // why: save bwdTC (i.e. objective functional gradient)
    IntermediateHistory,        // why: save for bwd
    RemainderHistory,           // why: save for bwd
    
    // Local scope, do this inline:
    EnergyEvolution,            // why: write energy to directory
    FourierSpectrumEvolution,   // why: write four to directory
    SolutionBranch,             // why: write branch to directory
    OptDiagnostics,             // why: write diag iterations to directory
    OptLineSearch               // why: write brent iterations to directory
    */  

    switch (storedDataType) {
        case InitialState: {
            std::ofstream file(paths.fOptimalInitialData, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case TerminalState: {
            std::ofstream file(paths.fOptimalTerminalData, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case BackwardInitialState: {
            std::ofstream file(paths.fBackwardSolution, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case IntermediateHistory:   
        case RemainderHistory: {
            std::ofstream file(appendTimeStep(paths.fForwardSolution, dCurrentT), std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
    }
}

void deleteData(const Pathnames &paths) { 
    //  delete after optimization due to temporary non-diagnostic use
    std::cout << "Deleting all temporary solution data...\n";
    std::filesystem::remove_all(paths.dirForwardSolution);
}