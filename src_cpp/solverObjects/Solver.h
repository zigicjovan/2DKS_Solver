#ifndef SOLVER_H
#define SOLVER_H

/*
#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "Timer.h"

#include <complex>
#include <vector>
#include <cstddef>

using Complex = std::complex<double>;

enum SolutionDataType {
    InitialState,               
    TerminalState,              
    BackwardInitialState,       // objective functional gradient
    IntermediateHistory,        
    RemainderHistory,           
};

class SolutionData {
private:
    std::vector<Complex> vForwardInitialState;
    std::vector<Complex> vForwardTerminalState;
    std::vector<Complex> vBackwardInitialState;
    std::vector<Complex> vBackwardTerminalState;

public:
    SolutionData(const Parameters& params, SolutionDataType storedDataType);
    
    std::size_t getSize() const;
    std::vector<Complex>& getData();
    Complex* getDataPointer(); // binary data                 
    std::vector<Complex>& setData(const std::vector<Complex>& stateData, int stateNumber);  
};
*/

#endif