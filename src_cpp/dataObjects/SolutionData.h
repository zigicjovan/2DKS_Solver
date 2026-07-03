#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"
#include "Pathnames.h"
#include "FFTWPlanner.h"

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
    Complex& operator[](std::size_t idx);
    Complex& atState(const Parameters& params, std::size_t p, int stateNumber);
    Complex* getStatePointer(const Parameters& params, int stateNumber);
    Complex* getDataPointer(); // binary data
    void setInitialEnergyL2(const Parameters& params);
    double getEnergyL2(const Parameters& params) const;
    double getEnergyH1(const Parameters& params) const;
    double getEnergyH2(const Parameters& params) const;
    std::vector<double> getRadialSpectrum(const Parameters& params) const;
};

#endif