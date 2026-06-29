#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"

#include <complex>
#include <vector>
#include <cstddef>

using Complex = std::complex<double>;

enum SolutionDataType {
    // Global scope:
    InitialState,               
    TerminalState,              
    BackwardInitialState,       // objective functional gradient
    IntermediateHistory,        
    RemainderHistory,           
    
    // // Local scope:
    // EnergyEvolution,            // write energy to directory
    // FourierSpectrumEvolution,   // write four to directory
    // SolutionBranch,             // load, manipulate, write branch to directory
    // OptDiagnostics,             // write diag iterations to directory
    // OptLineSearch               // write brent iterations to directory
};

class SolutionData {
private:
    std::vector<Complex> vData;

public:
    SolutionData(const Parameters& params, SolutionDataType storedDataType);
    
    std::size_t getSize() const;
    std::vector<Complex>& getData();
    Complex* getDataPointer(); // binary data                 
    std::vector<Complex>& setData(const std::vector<Complex>& stateData, int stateNumber);  
};

#endif