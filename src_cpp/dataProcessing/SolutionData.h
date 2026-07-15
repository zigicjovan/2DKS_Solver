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
    vector<complex<double>> _vData;
    const Parameters& _params;
    const Pathnames& _paths;

public:
    SolutionData(const Parameters& params, const Pathnames& paths, SolutionDataType storedDataType);
    
    size_t getSize() const;
    void getData(SolutionData&  stateData, size_t stateNumber = 0) const ;                
    void setData(const SolutionData&  stateData, size_t stateNumber = 0);  

    complex<double>& operator()(size_t i, size_t j, size_t stateNumber);
    complex<double>& operator[](size_t idx);
    complex<double>& atState(size_t p, size_t stateNumber);
    complex<double>* getStatePointer(size_t stateNumber);
    complex<double>* getDataPointer(); // for binary data
    void swapDataFrom(SolutionData& otherData) noexcept;
    void moveDataFrom(SolutionData& otherData) noexcept;

    void setInitialEnergyL2();
    double getEnergyL2() const;
    double getEnergyH1() const;
    double getEnergyH2() const;
    vector<double> getRadialSpectrum() const;

    filesystem::path appendTimeStep(const filesystem::path& path, double dCurrentT);
    void loadData(SolutionDataType storedDataType, double dCurrentT = 0.0);
    void saveData(SolutionDataType storedDataType, double dCurrentT = 0.0);
    void deleteData();
};

#endif