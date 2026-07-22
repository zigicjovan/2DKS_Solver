#ifndef SOLUTIONDATA_H
#define SOLUTIONDATA_H

#include "Parameters.h"
#include "Pathnames.h"
#include "FFTWPlanner.h"
#include "MPIContext.h"

#include <complex>
#include <vector>
#include <cstddef>
#include <filesystem>
#include <climits>
#include <mpi.h>

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
    const MPIContext& _mpi;

    void writeDistributedFile(const filesystem::path& filename) const;
    void readDistributedFile(const filesystem::path& filename);
    filesystem::path appendTimeStep(const filesystem::path& path, double dCurrentT);

public:
    SolutionData(const Parameters& params, const Pathnames& paths, const MPIContext& mpi, SolutionDataType storedDataType);
    
    void validateMPIConfiguration();
    size_t getSize() const;
    void getData(SolutionData&  stateData, size_t stateNumber = 0) const;                
    void setData(const SolutionData&  stateData, size_t stateNumber = 0);  

    SolutionData& operator=(const SolutionData& otherData);
    SolutionData& operator+=(const SolutionData& otherData);
    SolutionData& operator-=(const SolutionData& otherData);
    SolutionData& operator*=(double scalar);
    friend SolutionData operator+(SolutionData lhs, const SolutionData& rhs);
    friend SolutionData operator-(SolutionData lhs, const SolutionData& rhs);
    friend SolutionData operator*(SolutionData data, double scalar);
    friend SolutionData operator*(double scalar, SolutionData data);

    complex<double>& operator()(size_t i, size_t j, size_t stateNumber);
    complex<double>& operator[](size_t idx);
    const complex<double>& operator[](size_t idx) const;
    complex<double>& atState(size_t p, size_t stateNumber);
    complex<double>* getStatePointer(size_t stateNumber);
    const complex<double>* getStatePointer(size_t stateNumber) const;
    complex<double>* getDataPointer(); // for binary data
    void swapDataFrom(SolutionData& otherData) noexcept;
    void moveDataFrom(SolutionData& otherData) noexcept;

    void setInitialEnergyL2();
    double getInnerProductL2With(SolutionData& otherData) const;
    double getNormL2() const;
    double getEnergyL2() const;
    double getEnergyH1() const;
    double getEnergyH2() const;
    vector<double> getRadialSpectrum() const;

    void loadData(SolutionDataType storedDataType, double dCurrentT = 0.0);
    void saveData(SolutionDataType storedDataType, double dCurrentT = 0.0);
    void deleteData();
};

#endif