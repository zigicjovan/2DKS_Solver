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

SolutionData::SolutionData(const Parameters& params, const Pathnames& paths, SolutionDataType storedDataType) : _params(params), _paths(paths) {
    switch (storedDataType) {
        
        case InitialState: // fall through
        
        case TerminalState: // fall through
        
        case BackwardInitialState:
            _vData.resize(_params.getTotalGridSize() * 1);
            break;
    
        case IntermediateHistory:   
            _vData.resize(_params.getTotalGridSize() * _params.getNumericalStepsPerFile() );
            break;
        
        case RemainderHistory:
            _vData.resize(_params.getTotalGridSize() * ( _params.getNumericalSteps() % _params.getNumericalStepsPerFile()) );
            break;
    }
}

// Public functions
size_t SolutionData::getSize() const {
    return _vData.size();
}

void SolutionData::getData(SolutionData& stateData, size_t stateIndex) const {
    const size_t iOffset = stateIndex * _params.getTotalGridSize();
    copy(_vData.begin() + iOffset, _vData.begin() + iOffset + _params.getTotalGridSize(), stateData._vData.begin());
}

void SolutionData::setData(const SolutionData& stateData, size_t stateIndex) {
    const size_t iOffset = stateIndex * _params.getTotalGridSize();
    copy(stateData._vData.begin(), stateData._vData.begin() + _params.getTotalGridSize(), _vData.begin() + iOffset);
}

SolutionData& SolutionData::operator=(const SolutionData& otherData) {
    if (this != &otherData) {
        _vData = otherData._vData;
    }
    return *this;
}

SolutionData& SolutionData::operator+=(const SolutionData& otherData) {
    for (size_t i = 0; i < _vData.size(); ++i) {
        _vData[i] += otherData._vData[i];
    }
    return *this;
}

SolutionData& SolutionData::operator-=(const SolutionData& otherData) {
    for (size_t i = 0; i < _vData.size(); ++i) {
        _vData[i] -= otherData._vData[i];
    }
    return *this;
}

SolutionData& SolutionData::operator*=(double scalar) {
    for (complex<double>& value : _vData) {
        value *= scalar;
    }
    return *this;
}

SolutionData operator+(SolutionData lhs, const SolutionData& rhs) {
    lhs += rhs;
    return lhs;
}

SolutionData operator-(SolutionData lhs, const SolutionData& rhs) {
    lhs -= rhs;
    return lhs;
}

SolutionData operator*(SolutionData data, double scalar) {
    data *= scalar;
    return data;
}

SolutionData operator*(double scalar, SolutionData data) {
    data *= scalar;
    return data;
}

complex<double>& SolutionData::operator()(size_t i, size_t j, size_t stateNumber) {
    return _vData[stateNumber * _params.getTotalGridSize() + j * _params.getGridSize1() + i];
}

complex<double>& SolutionData::operator[](size_t idx) {
    return _vData[idx];
}

complex<double>& SolutionData::atState(size_t p, size_t stateNumber) {
    return _vData[stateNumber * _params.getTotalGridSize() + p];
}

complex<double>* SolutionData::getStatePointer(size_t stateNumber) {
    return _vData.data() + stateNumber * _params.getTotalGridSize();
}

complex<double>* SolutionData::getDataPointer() {
    return _vData.data();
}

void SolutionData::swapDataFrom(SolutionData& otherData) noexcept {
    _vData.swap(otherData._vData);
}

void SolutionData::moveDataFrom(SolutionData& otherData) noexcept {
    _vData = move(otherData._vData);
}

void SolutionData::setInitialEnergyL2() {
    double energyL2 = getEnergyL2();
    double normL2 = sqrt(energyL2);
    double scale = sqrt(_params.getInitialEnergy()) / normL2;
    size_t gridSize = _params.getTotalGridSize();

    for (size_t i = 0; i < gridSize; ++i) {
        _vData[i] *= scale;
    }
}

double SolutionData::getInnerProductL2With(SolutionData& otherData) const {
    double innerProductL2 = 0.0;
    size_t gridSize = _params.getTotalGridSize();

    for (size_t i = 0; i < gridSize; ++i) {
        innerProductL2 += real(conj(_vData[i]) * otherData._vData[i]);
    }
    return innerProductL2 * _params.getEnergyFactor();
}

double SolutionData::getNormL2() const {
    return sqrt(getEnergyL2());
}

double SolutionData::getEnergyL2() const {
    double energyL2 = 0.0;
    size_t gridSize = _params.getTotalGridSize();

    for (size_t i = 0; i < gridSize; ++i) {
        energyL2 += norm(_vData[i]);
    }
    return energyL2 * _params.getEnergyFactor();
}

double SolutionData::getEnergyH1() const {
    double energyH1 = 0.0;
    size_t gridSize = _params.getTotalGridSize();
    vector<double> vH1Weight = _params.getH1Weight();

    for (size_t i = 0; i < gridSize; ++i) {
        energyH1 += vH1Weight[i] * norm(_vData[i]);
    }
    return energyH1 * _params.getEnergyFactor();
}

double SolutionData::getEnergyH2() const {
    double energyH2 = 0.0;
    size_t gridSize = _params.getTotalGridSize();
    vector<double> vH2Weight = _params.getH2Weight();

    for (size_t i = 0; i < gridSize; ++i) {
        energyH2 += vH2Weight[i] * norm(_vData[i]);
    }
    return energyH2 * _params.getEnergyFactor();
}

vector<double> SolutionData::getRadialSpectrum() const {
    vector<double> spectrum(_params.getMaxRadialBin() + 1, 0.0);
    size_t gridSize = _params.getTotalGridSize();
    vector<size_t> radialBin = _params.getRadialBin();

    for (size_t i = 0; i < gridSize; ++i) {
        spectrum[radialBin[i]] += norm(_vData[i]);
    }

    for (double& value : spectrum) {
        value *= 0.5;
    }

    return spectrum;
}

void SolutionData::loadData(SolutionDataType storedDataType, double dCurrentT) {   
    switch (storedDataType) {

        case InitialState: {
            ifstream file(_paths.getInitialDataFile(), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }

        case TerminalState: {
            ifstream file(_paths.getTerminalDataFile(), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }

        case BackwardInitialState: {
            ifstream file(_paths.getBackwardSolutionFile(), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case IntermediateHistory: // fall through  
        
        case RemainderHistory: {
            ifstream file(appendTimeStep(_paths.getForwardSolutionFile(), dCurrentT), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));
            break;
        }
    }
}

void SolutionData::saveData(SolutionDataType storedDataType, double dCurrentT) {   
    switch (storedDataType) {

        case InitialState: {
            ofstream file(_paths.getInitialDataFile(), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case TerminalState: {
            ofstream file(_paths.getTerminalDataFile(), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case BackwardInitialState: {
            ofstream file(_paths.getBackwardSolutionFile(), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }

        case IntermediateHistory: // fall through 

        case RemainderHistory: {
            ofstream file(appendTimeStep(_paths.getForwardSolutionFile(), dCurrentT), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); 
            break;
        }
    }
}

void SolutionData::deleteData() { 
    filesystem::remove_all(_paths.getDirForwardSolution());
    filesystem::create_directories(_paths.getDirForwardSolution());
}