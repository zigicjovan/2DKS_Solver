#include "SolutionData.h"
#include "Parameters.h"

#include <complex>
#include <vector>

#include <iostream>  

using Complex = std::complex<double>;

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
            std::cout << "Intermediate States: " << vData.size() / params.iTotalGridSize << "\n";
            break;
        case RemainderHistory:
            vData.resize(params.iTotalGridSize * ( params.iGetNumericalSteps() % params.iGetNumericalStepsPerFile()) );
            std::cout << "Remainder States: " << vData.size() / params.iTotalGridSize << "\n";
            break;
    }
}

// Public functions

std::size_t SolutionData::getSize() const {
    return vData.size();
}

std::vector<Complex>& SolutionData::getData() {
    return vData;
}

void SolutionData::setData(const std::vector<Complex>& stateData, int stateIndex) {
    // stateData replaces part of vData at stateIndex
    std::copy(stateData.begin(), stateData.end(), vData.begin() + stateIndex);
}

Complex& SolutionData::operator()(const Parameters& params, std::size_t i, std::size_t j, int stateNumber) {
    return vData[stateNumber * params.iTotalGridSize + j * params.iGridSize1 + i];
}

Complex& SolutionData::operator[](std::size_t idx) {
    return vData[idx];
}

Complex& SolutionData::atState(const Parameters& params, std::size_t p, int stateNumber) {
    return vData[stateNumber * params.iTotalGridSize + p];
}

Complex* SolutionData::getStatePointer(const Parameters& params, int stateNumber) {
    return vData.data() + stateNumber * params.iTotalGridSize;
}

Complex* SolutionData::getDataPointer() {
    return vData.data();
}

void SolutionData::setInitialEnergyL2(const Parameters& params) {
    double energyL2 = getEnergyL2(params);
    double normL2 = std::sqrt(energyL2);
    double scale = std::sqrt(params.dInitialEnergy) / normL2;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        vData[p] *= scale;
}

double SolutionData::getEnergyL2(const Parameters& params) const {
    double energyL2 = 0.0;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        energyL2 += std::norm(vData[p]);
    return energyL2 * params.dEnergyFactor;
}

double SolutionData::getEnergyH1(const Parameters& params) const {
    double energyH1 = 0.0;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        energyH1 += params.vH1Weight[p] * std::norm(vData[p]);
    return energyH1 * params.dEnergyFactor;
}

double SolutionData::getEnergyH2(const Parameters& params) const {
    double energyH2 = 0.0;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        energyH2 += params.vH2Weight[p] * std::norm(vData[p]);
    return energyH2 * params.dEnergyFactor;
}

std::vector<double> SolutionData::getRadialSpectrum(const Parameters& params) const {
    std::vector<double> spectrum(params.iMaxRadialBin + 1, 0.0);
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        spectrum[params.vRadialBin[p]] += std::norm(vData[p]);
    for (double& value : spectrum)
        value *= 0.5;
    return spectrum;
}