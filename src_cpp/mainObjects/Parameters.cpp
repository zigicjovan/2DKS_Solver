#include "Parameters.h"
#include <cmath> 
#include <string>

#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <array>

using Complex = std::complex<double>;
std::size_t getIndex(std::size_t i, std::size_t j, std::size_t N) {
    return i * N + j;
}

// Private functions

std::vector<double> Parameters::setPhysicalSpace1() {
    std::vector<double> gridpoints(iGridSize1);
    for (std::size_t j = 0; j < iGridSize1; ++j)
        gridpoints[j] = dDomainSize1 * static_cast<double>(j) / static_cast<double>(iGridSize1);
    return gridpoints;
}

std::vector<double> Parameters::setPhysicalSpace2() {
    std::vector<double> gridpoints(iGridSize2);
    for (std::size_t j = 0; j < iGridSize2; ++j)
        gridpoints[j] = dDomainSize2 * static_cast<double>(j) / static_cast<double>(iGridSize2);
    return gridpoints;
}

std::vector<double> Parameters::setFourierModesNonlinear1() {
    std::vector<double> wavenumbers(iGridSize1);
    for (std::size_t j = 0; j < iGridSize1 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j);
    wavenumbers[iGridSize1 / 2] = 0.0;
    for (std::size_t j = iGridSize1 / 2 + 1; j < iGridSize1; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * (static_cast<double>(j) - static_cast<double>(iGridSize1));
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesNonlinear2() {
    std::vector<double> wavenumbers(iGridSize2);
    for (std::size_t j = 0; j < iGridSize2 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j);
    wavenumbers[iGridSize2 / 2] = 0.0;
    for (std::size_t j = iGridSize2 / 2 + 1; j < iGridSize2; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * (static_cast<double>(j) - static_cast<double>(iGridSize2));
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesLinear1() {
    std::vector<double> wavenumbers(iGridSize1);
    for (std::size_t j = 0; j <= iGridSize1 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j);

    for (std::size_t j = iGridSize1 / 2 + 1; j < iGridSize1; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * (static_cast<double>(j) - static_cast<double>(iGridSize1));
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesLinear2() {
    std::vector<double> wavenumbers(iGridSize2);
    for (std::size_t j = 0; j <= iGridSize2 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j);

    for (std::size_t j = iGridSize2 / 2 + 1; j < iGridSize2; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * (static_cast<double>(j) - static_cast<double>(iGridSize2));
    return wavenumbers;
}

// Public functions

std::size_t Parameters::iGetNumericalSteps() const {   
    return static_cast<std::size_t>(std::round(dTimeWindow / dTimeStep));
}

std::size_t Parameters::iGetNumericalStepsPerFile() const {   
    return static_cast<std::size_t>(std::round( 4e7 / (iTotalGridSize) ));
}

void Parameters::getMathematicalOperators() {
    dEnergyFactor = dDomainFactor1 * dDomainFactor2 * (2.0 * dPI) * (2.0 * dPI) / (static_cast<double>(iTotalGridSize) * iTotalGridSize);
    
    _vGridpoints1 = setPhysicalSpace1();
    _vGridpoints2 = setPhysicalSpace2();

    const double _dSpaceStep1 = dDomainSize1 / static_cast<double>(iGridSize1);
    const double _dSpaceStep2 = dDomainSize2 / static_cast<double>(iGridSize2);
    dSpaceStep = std::min(_dSpaceStep1, _dSpaceStep2);

    vWavenumbersNonlinear1 = setFourierModesNonlinear1();
    vWavenumbersNonlinear2 = setFourierModesNonlinear2();

    _vWavenumbersLinear1 = setFourierModesLinear1();
    _vWavenumbersLinear2 = setFourierModesLinear2();

    vGrid1.resize(iTotalGridSize);
    vGrid2.resize(iTotalGridSize);

    _vSpectrumNonlinear1.resize(iTotalGridSize);
    _vSpectrumNonlinear2.resize(iTotalGridSize);
    _vSpectrumLinear1.resize(iTotalGridSize);
    _vSpectrumLinear2.resize(iTotalGridSize);

    vH1Weight.resize(iTotalGridSize);
    vH2Weight.resize(iTotalGridSize);
    double maxH1 = 0.0;
    double maxH2 = 0.0;

    vLinearOperator.resize(iTotalGridSize);
    vLaplaceOperator.resize(iTotalGridSize);
    vDifferentialOperator1.resize(iTotalGridSize);
    vDifferentialOperator2.resize(iTotalGridSize);
    vDealiasingOperator.resize(iTotalGridSize);
    const double kCut1 = (2.0 / 3.0) * (static_cast<double>(iGridSize1) / 2.0);
    const double kCut2 = (2.0 / 3.0) * (static_cast<double>(iGridSize2) / 2.0);

    double maxRadius = 0.0;
    vRadialBin.resize(iTotalGridSize);

    for (std::size_t i = 0; i < iGridSize2; ++i) {
        for (std::size_t j = 0; j < iGridSize1; ++j) {
            const std::size_t p = getIndex(i, j, iGridSize1);
            vGrid1[p] = _vGridpoints1[j];
            vGrid2[p] = _vGridpoints2[i];
            _vSpectrumNonlinear1[p] = vWavenumbersNonlinear1[j];
            _vSpectrumNonlinear2[p] = vWavenumbersNonlinear2[i];
            _vSpectrumLinear1[p] = _vWavenumbersLinear1[j];
            _vSpectrumLinear2[p] = _vWavenumbersLinear2[i];
            
            const double _dSpectralRadiusSquared = _vSpectrumLinear1[p] * _vSpectrumLinear1[p] + _vSpectrumLinear2[p] * _vSpectrumLinear2[p];
            const double _dBilaplacianValue = _dSpectralRadiusSquared * _dSpectralRadiusSquared;
            vH1Weight[p] = 1.0 + _dSpectralRadiusSquared;
            vH2Weight[p] = vH1Weight[p] * vH1Weight[p];
            maxH1 = std::max(maxH1, vH1Weight[p]);
            maxH2 = std::max(maxH2, vH2Weight[p]);
            
            vLinearOperator[p] = -_dSpectralRadiusSquared + _dBilaplacianValue;
            vLaplaceOperator[p] = -_dSpectralRadiusSquared;
            vDifferentialOperator1[p] = Imaginary * _vSpectrumNonlinear1[p];
            vDifferentialOperator2[p] = Imaginary * _vSpectrumNonlinear2[p];

            const bool keepMode = std::abs(vWavenumbersNonlinear1[j]) <= kCut1 && std::abs(vWavenumbersNonlinear2[i]) <= kCut2;
            vDealiasingOperator[p] = keepMode ? Complex{1.0, 0.0} : Complex{0.0, 0.0};

            const double k1 = _vWavenumbersLinear1[j];
            const double k2 = _vWavenumbersLinear2[i];
            const double radius = std::sqrt(k1 * k1 + k2 * k2);
            vRadialBin[p] = static_cast<std::size_t>(std::round(radius));
            maxRadius = std::max(maxRadius, radius);
        }   
    }
    iMaxRadialBin = static_cast<std::size_t>(std::ceil(maxRadius));
}