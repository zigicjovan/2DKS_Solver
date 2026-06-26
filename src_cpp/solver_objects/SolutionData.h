#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"

#include <complex>
#include <vector>
#include <cstddef>

using Complex = std::complex<double>;

enum class SolutionDataType
{
    // Global scope:
    InitialState,               // why: for storing optIC
    TerminalState,              // why: for storing optTC
    BackwardInitialState,       // why: for storing bwdTC (i.e. objective functional gradient)
    IntermediateHistory,        // why: for storing fwd
    RemainderHistory,           // why: for storing fwd
    
    // // Local scope:
    // EnergyEvolution,            // why: write energy to directory
    // FourierSpectrumEvolution,   // why: write four to directory
    // SolutionBranch,             // why: load, manipulate, write branch to directory
    // OptDiagnostics,             // why: write diag iterations to directory
    // OptLineSearch               // why: write brent iterations to directory
};

class SolutionData
{
private:
    std::vector<Complex> vec_Data;

public:
    SolutionData(const Parameters &params, SolutionDataType storedDataType);
    
    std::size_t getSize() const;
    std::vector<Complex> &getData();
    Complex *getDataPointer();                  // why: required pointers (not reference) to write binary data to directory      
    // const Complex *getDataPointer() const;   // why: required pointer (not reference) to write binary data to directory 

    // std::vector<Complex> &setData(const std::vector<Complex> &state);  
};

#endif