#include "Pathnames.h"
#include "SolutionData.h"

#include <filesystem>
#include <fstream>
#include <iostream> 

std::filesystem::path appendTimeStep(const std::filesystem::path& path, int currentT)
{
    return path.parent_path() /
           (path.stem().string()
            + "_" + std::to_string(currentT)
            + path.extension().string());
}

void loadData(const Pathnames &paths, SolutionData &storedData, SolutionDataType storedDataType, int currentT)
{   /*
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

    switch (storedDataType)
    {
        case SolutionDataType::InitialState:
        {
            std::ifstream file(paths.file_OptimalInitialData, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex));// why: save as binary
            break;
        }
        case SolutionDataType::TerminalState:
        {
            std::ifstream file(paths.file_OptimalTerminalData, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case SolutionDataType::BackwardInitialState:
        {
            std::ifstream file(paths.file_BackwardSolution, std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case SolutionDataType::IntermediateHistory:   
        case SolutionDataType::RemainderHistory:
        {
            std::ifstream file(appendTimeStep(paths.file_ForwardSolution, currentT), std::ios::binary);
            file.read(reinterpret_cast<char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
    }
}

void saveData(const Pathnames &paths, SolutionData &storedData, SolutionDataType storedDataType, int currentT)
{   /*
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

    switch (storedDataType)
    {
        case SolutionDataType::InitialState:
        {
            std::ofstream file(paths.file_OptimalInitialData, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case SolutionDataType::TerminalState:
        {
            std::ofstream file(paths.file_OptimalTerminalData, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case SolutionDataType::BackwardInitialState:
        {
            std::ofstream file(paths.file_BackwardSolution, std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
        case SolutionDataType::IntermediateHistory:   
        case SolutionDataType::RemainderHistory:
        {
            std::ofstream file(appendTimeStep(paths.file_ForwardSolution, currentT), std::ios::binary);
            file.write(reinterpret_cast<const char*>(storedData.getDataPointer()), storedData.getSize() * sizeof(Complex)); // why: save as binary
            break;
        }
    }
}

void deleteData(const Pathnames &paths)
{   /*
    SolutionDataTypes:
    // Global scope:
    IntermediateHistory,        // why: delete after optimization due to temporary non-diagnostic use
    RemainderHistory,           // why: delete after optimization due to temporary non-diagnostic use
    */ 
    std::cout << "Deleting all temporary solution data...\n";
    std::filesystem::remove_all(paths.dir_ForwardSolution);
}