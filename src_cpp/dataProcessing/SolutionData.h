#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"
#include "Pathnames.h"
#include "FFTWPlanner.h"

#include <complex>
#include <vector>
#include <cstddef>
#include <filesystem>

using namespace std;

enum SolutionDataType {
    InitialState,               
    TerminalState,              
    BackwardInitialState,       
    IntermediateHistory,        
    RemainderHistory,           
};

class SolutionData {
private:
    vector<complex<double>> vData;

public:
    SolutionData(const Parameters& params, SolutionDataType storedDataType);
    size_t getSize() const;
    vector<complex<double>>& getData();                 
    void setData(const vector<complex<double>>& stateData, int stateNumber);  

    complex<double>& operator()(const Parameters& params, size_t i, size_t j, int stateNumber);
    complex<double>& operator[](size_t idx);
    complex<double>& atState(const Parameters& params, size_t p, int stateNumber);
    complex<double>* getStatePointer(const Parameters& params, int stateNumber);
    complex<double>* getDataPointer(); // binary data
    void setInitialEnergyL2(const Parameters& params);
    double getEnergyL2(const Parameters& params) const;
    double getEnergyH1(const Parameters& params) const;
    double getEnergyH2(const Parameters& params) const;
    vector<double> getRadialSpectrum(const Parameters& params) const;

    filesystem::path appendTimeStep(const filesystem::path& path, double dCurrentT);
    void loadData(const Pathnames& paths, SolutionDataType storedDataType, double dCurrentT = 0.0);
    void saveData(const Pathnames& paths, SolutionDataType storedDataType, double dCurrentT = 0.0);
    void deleteData(const Pathnames& paths);
};

#endif