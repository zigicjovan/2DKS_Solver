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
SolutionData::SolutionData(const Parameters& params, SolutionDataType storedDataType) {
    switch (storedDataType) {
        case InitialState:
        case TerminalState:
        case BackwardInitialState:
            vData.resize(params.iTotalGridSize * 1);
            break;
        case IntermediateHistory:   
            vData.resize(params.iTotalGridSize * params.iGetNumericalStepsPerFile() );
            cout << "Intermediate States: " << vData.size() / params.iTotalGridSize << "\n";
            break;
        case RemainderHistory:
            vData.resize(params.iTotalGridSize * ( params.iGetNumericalSteps() % params.iGetNumericalStepsPerFile()) );
            cout << "Remainder States: " << vData.size() / params.iTotalGridSize << "\n";
            break;
    }
}

// Public functions
size_t SolutionData::getSize() const {
    return vData.size();
}

vector<complex<double>>& SolutionData::getData() {
    return vData;
}

void SolutionData::setData(const vector<complex<double>>& stateData, int stateIndex) {
    
    // stateData replaces part of vData at stateIndex
    copy(stateData.begin(), stateData.end(), vData.begin() + stateIndex);
}

complex<double>& SolutionData::operator()(const Parameters& params, size_t i, size_t j, int stateNumber) {
    return vData[stateNumber * params.iTotalGridSize + j * params.iGridSize1 + i];
}

complex<double>& SolutionData::operator[](size_t idx) {
    return vData[idx];
}

complex<double>& SolutionData::atState(const Parameters& params, size_t p, int stateNumber) {
    return vData[stateNumber * params.iTotalGridSize + p];
}

complex<double>* SolutionData::getStatePointer(const Parameters& params, int stateNumber) {
    return vData.data() + stateNumber * params.iTotalGridSize;
}

complex<double>* SolutionData::getDataPointer() {
    return vData.data();
}

void SolutionData::setInitialEnergyL2(const Parameters& params) {
    double energyL2 = getEnergyL2(params);
    double normL2 = sqrt(energyL2);
    double scale = sqrt(params.dInitialEnergy) / normL2;
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        vData[i] *= scale;
}

double SolutionData::getEnergyL2(const Parameters& params) const {
    double energyL2 = 0.0;
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        energyL2 += norm(vData[i]);
    return energyL2 * params.dEnergyFactor;
}

double SolutionData::getEnergyH1(const Parameters& params) const {
    double energyH1 = 0.0;
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        energyH1 += params.vH1Weight[i] * norm(vData[i]);
    return energyH1 * params.dEnergyFactor;
}

double SolutionData::getEnergyH2(const Parameters& params) const {
    double energyH2 = 0.0;
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        energyH2 += params.vH2Weight[i] * norm(vData[i]);
    return energyH2 * params.dEnergyFactor;
}

vector<double> SolutionData::getRadialSpectrum(const Parameters& params) const {
    vector<double> spectrum(params.iMaxRadialBin + 1, 0.0);
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        spectrum[params.vRadialBin[i]] += norm(vData[i]);
    for (double& value : spectrum)
        value *= 0.5;
    return spectrum;
}

void SolutionData::loadData(const Pathnames& paths, SolutionDataType storedDataType, double dCurrentT) {   
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
            ifstream file(paths.fOptimalInitialData, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>));// why: save as binary
            break;
        }

        case TerminalState: {
            ifstream file(paths.fOptimalTerminalData, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }

        case BackwardInitialState: {
            ifstream file(paths.fBackwardSolution, ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }

        case IntermediateHistory: // fall through  
        
        case RemainderHistory: {
            ifstream file(appendTimeStep(paths.fForwardSolution, dCurrentT), ios::binary);
            file.read(reinterpret_cast<char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }
    }
}

void SolutionData::saveData(const Pathnames& paths, SolutionDataType storedDataType, double dCurrentT) {   
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
            ofstream file(paths.fOptimalInitialData, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }

        case TerminalState: {
            ofstream file(paths.fOptimalTerminalData, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }

        case BackwardInitialState: {
            ofstream file(paths.fBackwardSolution, ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }

        case IntermediateHistory: // fall through 

        case RemainderHistory: {
            ofstream file(appendTimeStep(paths.fForwardSolution, dCurrentT), ios::binary);
            file.write(reinterpret_cast<const char*>(this->getDataPointer()), this->getSize() * sizeof(complex<double>)); // why: save as binary
            break;
        }
    }
}

void SolutionData::deleteData(const Pathnames &paths) { 

    //  delete after optimization due to temporary non-diagnostic use
    cout << "Deleting all temporary solution data...\n";
    filesystem::remove_all(paths.dirForwardSolution);
}