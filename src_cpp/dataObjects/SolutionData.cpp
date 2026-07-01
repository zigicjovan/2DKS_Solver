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

Complex* SolutionData::getStatePointer(const Parameters& params, int stateNumber) {
    return vData.data() + stateNumber * params.iTotalGridSize;
}

Complex* SolutionData::getDataPointer() {
    return vData.data();
}

void SolutionData::setInitialEnergyL2(const Parameters& params) {
    double energyL2 = 0.0;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        energyL2 += std::norm(vData[p]);
    energyL2 *= params.dDomainFactor1 * params.dDomainFactor2
            * std::pow(2.0 * dPI, 2) / std::pow(static_cast<double>(params.iTotalGridSize), 2);
    
    double normL2 = std::sqrt(energyL2);
    double scale = std::sqrt(params.dInitialEnergy) / normL2;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        vData[p] *= scale;
}

double SolutionData::getEnergyL2(const Parameters& params, const Complex* vFourierState) {
    double energy = 0.0;
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        energy += std::norm(vFourierState[p]);
    energy *= params.dDomainFactor1 * params.dDomainFactor2
            * std::pow(2.0 * dPI, 2) / std::pow(static_cast<double>(params.iTotalGridSize), 2);
    return energy;
}