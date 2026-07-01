#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"

#include <complex>
#include <vector>
#include <cstddef>

using Complex = std::complex<double>;

enum SolutionDataType {
    InitialState,               
    TerminalState,              
    BackwardInitialState,       
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
    void setData(const std::vector<Complex>& stateData, int stateNumber);  
    Complex& operator()(const Parameters& params, std::size_t i, std::size_t j, int stateNumber);
    Complex* getStatePointer(const Parameters& params, int stateNumber);
    Complex* getDataPointer(); // binary data
    void setInitialEnergyL2(const Parameters& params);
    double getEnergyL2(const Parameters& params, const Complex* vState);
};

#endif