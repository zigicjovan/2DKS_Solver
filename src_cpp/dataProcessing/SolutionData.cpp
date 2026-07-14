#include "SolutionData.h"
#include "Parameters.h"

#include <complex>
#include <vector> 
#include <filesystem>
#include <fstream>
#include <iostream> 

using namespace std;

filesystem::path SolutionData::appendTimeStep(const filesystem::path& path, double dCurrentT) {
    return path.parent_path() / (path.stem().string() + "_" + to_string(dCurrentT) + path.extension().string());
}

// Fourier space solution data
SolutionData::SolutionData(const Parameters& params, const Pathnames& paths, SolutionDataType storedDataType) : _params(params), _paths(paths) {
    switch (storedDataType) {
        case InitialState:
        case TerminalState:
        case BackwardInitialState:
            _vData.resize(_params.iTotalGridSize * 1);
            break;
        case IntermediateHistory:   
            _vData.resize(_params.iTotalGridSize * _params.iGetNumericalStepsPerFile() );
            cout << "Intermediate States: " << _vData.size() / _params.iTotalGridSize << "\n";
            break;
        case RemainderHistory:
            _vData.resize(_params.iTotalGridSize * ( _params.iGetNumericalSteps() % _params.iGetNumericalStepsPerFile()) );
            cout << "Remainder States: " << _vData.size() / _params.iTotalGridSize << "\n";
            break;
    }
}

// Public functions
size_t SolutionData::getSize() const {
    return _vData.size();
}

vector<complex<double>>& SolutionData::getData() {
    return _vData;
}

void SolutionData::setData(const SolutionData& stateData, size_t stateIndex) {
    const size_t iOffset = stateIndex * _params.iTotalGridSize;
    copy(stateData._vData.begin(), stateData._vData.end(), _vData.begin() + iOffset);
}

complex<double>& SolutionData::operator()(size_t i, size_t j, size_t stateNumber) {
    return _vData[stateNumber * _params.iTotalGridSize + j * _params.iGridSize1 + i];
}

complex<double>& SolutionData::operator[](size_t idx) {
    return _vData[idx];
}

complex<double>& SolutionData::atState(size_t p, size_t stateNumber) {
    return _vData[stateNumber * _params.iTotalGridSize + p];
}

complex<double>* SolutionData::getStatePointer(size_t stateNumber) {
    return _vData.data() + stateNumber * _params.iTotalGridSize;
}

complex<double>* SolutionData::getDataPointer() {
    return _vData.data();
}

void SolutionData::swapDataFrom(SolutionData& otherData) noexcept
{
    _vData.swap(otherData._vData);
}

void SolutionData::moveDataFrom(SolutionData& otherData) noexcept
{
    _vData = move(otherData._vData);
}

void SolutionData::setInitialEnergyL2() {
    double energyL2 = getEnergyL2();
    double normL2 = sqrt(energyL2);
    double scale = sqrt(_params.dInitialEnergy) / normL2;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        _vData[i] *= scale;
}

double SolutionData::getEnergyL2() const {
    double energyL2 = 0.0;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        energyL2 += norm(_vData[i]);
    return energyL2 * _params.dEnergyFactor;
}

double SolutionData::getEnergyH1() const {
    double energyH1 = 0.0;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        energyH1 += _params.vH1Weight[i] * norm(_vData[i]);
    return energyH1 * _params.dEnergyFactor;
}

double SolutionData::getEnergyH2() const {
    double energyH2 = 0.0;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        energyH2 += _params.vH2Weight[i] * norm(_vData[i]);
    return energyH2 * _params.dEnergyFactor;
}

vector<double> SolutionData::getRadialSpectrum() const {
    vector<double> spectrum(_params.iMaxRadialBin + 1, 0.0);
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        spectrum[_params.vRadialBin[i]] += norm(_vData[i]);
    for (double& value : spectrum)
        value *= 0.5;
    return spectrum;
}

void SolutionData::loadData(SolutionDataType storedDataType, double dCurrentT) {   
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
            ifstream file(_paths.fOptimalInitialData, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }

        case TerminalState: {
            ifstream file(_paths.fOptimalTerminalData, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }

        case BackwardInitialState: {
            ifstream file(_paths.fBackwardSolution, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case IntermediateHistory: // fall through  
        
        case RemainderHistory: {
            ifstream file(appendTimeStep(_paths.fForwardSolution, dCurrentT), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }
    }
}

void SolutionData::saveData(SolutionDataType storedDataType, double dCurrentT) {   
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
            ofstream file(_paths.fOptimalInitialData, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case TerminalState: {
            ofstream file(_paths.fOptimalTerminalData, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case BackwardInitialState: {
            ofstream file(_paths.fBackwardSolution, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case IntermediateHistory: // fall through 

        case RemainderHistory: {
            ofstream file(appendTimeStep(_paths.fForwardSolution, dCurrentT), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }
    }
}

void SolutionData::deleteData() { 

    //  delete after optimization due to temporary non-diagnostic use
    cout << "Deleting all temporary solution data...\n";
    filesystem::remove_all(_paths.dirForwardSolution);
}