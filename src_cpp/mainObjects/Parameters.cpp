#include "Parameters.h"
#include <cmath> 
#include <string>

#include <vector>
#include <complex>
#include <algorithm>
#include <array>

using Complex = std::complex<double>;
int getIndex(int i, int j, int N) {
    return i * N + j;
}

// Private functions

std::vector<double> Parameters::setPhysicalSpace1() {
    std::vector<double> gridpoints(iGridSize1);
    for (int j = 0; j < iGridSize1; ++j)
        gridpoints[j] = dDomainSize1 * static_cast<double>(j) / static_cast<double>(iGridSize1);
    return gridpoints;
}

std::vector<double> Parameters::setPhysicalSpace2() {
    std::vector<double> gridpoints(iGridSize2);
    for (int j = 0; j < iGridSize2; ++j)
        gridpoints[j] = dDomainSize2 * static_cast<double>(j) / static_cast<double>(iGridSize2);
    return gridpoints;
}

std::vector<double> Parameters::setFourierModesNonlinear1() {
    std::vector<double> wavenumbers(iGridSize1);
    for (int j = 0; j < iGridSize1 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j);
    wavenumbers[iGridSize1 / 2] = 0.0;
    for (int j = iGridSize1 / 2 + 1; j < iGridSize1; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j - iGridSize1);
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesNonlinear2() {
    std::vector<double> wavenumbers(iGridSize2);
    for (int j = 0; j < iGridSize2 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j);
    wavenumbers[iGridSize2 / 2] = 0.0;
    for (int j = iGridSize2 / 2 + 1; j < iGridSize2; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j - iGridSize2);
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesLinear1() {
    std::vector<double> wavenumbers(iGridSize1);
    for (int j = 0; j <= iGridSize1 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j);
    wavenumbers[iGridSize1 / 2] = 0.0;
    for (int j = iGridSize1 / 2 + 1; j < iGridSize1; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize1 * static_cast<double>(j - iGridSize1);
    return wavenumbers;
}

std::vector<double> Parameters::setFourierModesLinear2() {
    std::vector<double> wavenumbers(iGridSize2);
    for (int j = 0; j <= iGridSize2 / 2; ++j) 
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j);
    wavenumbers[iGridSize2 / 2] = 0.0;
    for (int j = iGridSize2 / 2 + 1; j < iGridSize2; ++j)
        wavenumbers[j] = 2.0 * dPI / dDomainSize2 * static_cast<double>(j - iGridSize2);
    return wavenumbers;
}

// Public functions

int Parameters::iGetNumericalSteps() const {   
    return static_cast<int>(std::round(dTimeWindow / dTimeStep));
}

int Parameters::iGetNumericalStepsPerFile() const {   
    return static_cast<int>(std::round( 4e7 / (iGridSize1 * iGridSize2) ));
}

void Parameters::getMathematicalOperators() {
    _vGridpoints1 = setPhysicalSpace1();
    _vGridpoints2 = setPhysicalSpace2();

    const double _dSpaceStep1 = dDomainSize1 / static_cast<double>(iGridSize1);
    const double _dSpaceStep2 = dDomainSize2 / static_cast<double>(iGridSize2);
    dSpaceStep = std::min(_dSpaceStep1, _dSpaceStep2);

    _vWavenumbersNonlinear1 = setFourierModesNonlinear1();
    _vWavenumbersNonlinear2 = setFourierModesNonlinear2();

    _vWavenumbersLinear1 = setFourierModesLinear1();
    _vWavenumbersLinear2 = setFourierModesLinear2();

    _vGrid1.resize(iGridSize1 * iGridSize2);
    _vGrid2.resize(iGridSize1 * iGridSize2);

    _vSpectrumNonlinear1.resize(iGridSize1 * iGridSize2);
    _vSpectrumNonlinear2.resize(iGridSize1 * iGridSize2);
    _vSpectrumLinear1.resize(iGridSize1 * iGridSize2);
    _vSpectrumLinear2.resize(iGridSize1 * iGridSize2);

    for (int i = 0; i < iGridSize2; ++i) {
        for (int j = 0; j < iGridSize1; ++j) {
            const int p = getIndex(i, j, iGridSize1);
            _vGrid1[p] = _vGridpoints1[j];
            _vGrid2[p] = _vGridpoints2[i];
            _vSpectrumNonlinear1[p] = _vWavenumbersNonlinear1[j];
            _vSpectrumNonlinear2[p] = _vWavenumbersNonlinear2[i];
            _vSpectrumLinear1[p] = _vWavenumbersLinear1[j];
            _vSpectrumLinear2[p] = _vWavenumbersLinear2[i];
            const double _dLaplacianValue = _vSpectrumLinear1[p] * _vSpectrumLinear1[p] + _vSpectrumLinear2[p] * _vSpectrumLinear2[p];
            const double _dBilaplacianValue = _dLaplacianValue * _dLaplacianValue;
            vLinearOperator[p] = -_dLaplacianValue + _dBilaplacianValue;
            vLaplaceOperator[p] = -_dLaplacianValue;
            vDifferentialOperator1[p] = Imaginary * _vSpectrumNonlinear1[p];
            vDifferentialOperator2[p] = Imaginary * _vSpectrumNonlinear2[p];
        }   
    }
}