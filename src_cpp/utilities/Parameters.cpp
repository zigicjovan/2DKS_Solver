#include "Parameters.h"
#include <cmath> 
#include <string>

#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <array>

using namespace std;

size_t getIndex(size_t i, size_t j, size_t N) {
    return i * N + j;
}

// Private functions
vector<double> Parameters::setPhysicalSpace1() {
    vector<double> gridpoints(iGridSize1);
    for (size_t i = 0; i < iGridSize1; ++i)
        gridpoints[i] = dDomainSize1 * static_cast<double>(i) / static_cast<double>(iGridSize1);
    return gridpoints;
}

vector<double> Parameters::setPhysicalSpace2() {
    vector<double> gridpoints(iGridSize2);
    for (size_t i = 0; i < iGridSize2; ++i)
        gridpoints[i] = dDomainSize2 * static_cast<double>(i) / static_cast<double>(iGridSize2);
    return gridpoints;
}

vector<double> Parameters::setFourierModesNonlinear1() {
    vector<double> wavenumbers(iGridSize1);
    for (size_t i = 0; i < iGridSize1 / 2; ++i) 
        wavenumbers[i] = 2.0 * dPI / dDomainSize1 * static_cast<double>(i);
    wavenumbers[iGridSize1 / 2] = 0.0;
    for (size_t i = iGridSize1 / 2 + 1; i < iGridSize1; ++i)
        wavenumbers[i] = 2.0 * dPI / dDomainSize1 * (static_cast<double>(i) - static_cast<double>(iGridSize1));
    return wavenumbers;
}

vector<double> Parameters::setFourierModesNonlinear2() {
    vector<double> wavenumbers(iGridSize2);
    for (size_t i = 0; i < iGridSize2 / 2; ++i) 
        wavenumbers[i] = 2.0 * dPI / dDomainSize2 * static_cast<double>(i);
    wavenumbers[iGridSize2 / 2] = 0.0;
    for (size_t i = iGridSize2 / 2 + 1; i < iGridSize2; ++i)
        wavenumbers[i] = 2.0 * dPI / dDomainSize2 * (static_cast<double>(i) - static_cast<double>(iGridSize2));
    return wavenumbers;
}

vector<double> Parameters::setFourierModesLinear1() {
    vector<double> wavenumbers(iGridSize1);
    for (size_t i = 0; i <= iGridSize1 / 2; ++i) 
        wavenumbers[i] = 2.0 * dPI / dDomainSize1 * static_cast<double>(i);

    for (size_t i = iGridSize1 / 2 + 1; i < iGridSize1; ++i)
        wavenumbers[i] = 2.0 * dPI / dDomainSize1 * (static_cast<double>(i) - static_cast<double>(iGridSize1));
    return wavenumbers;
}

vector<double> Parameters::setFourierModesLinear2() {
    vector<double> wavenumbers(iGridSize2);
    for (size_t i = 0; i <= iGridSize2 / 2; ++i) 
        wavenumbers[i] = 2.0 * dPI / dDomainSize2 * static_cast<double>(i);

    for (size_t i = iGridSize2 / 2 + 1; i < iGridSize2; ++i)
        wavenumbers[i] = 2.0 * dPI / dDomainSize2 * (static_cast<double>(i) - static_cast<double>(iGridSize2));
    return wavenumbers;
}

// Public functions
Parameters::Parameters(int argc, char* argv[]) {

    strInitialGuessName = argv[1];
    iGridSize1 = stoi(argv[2]);
    iGridSize2 = stoi(argv[3]);
    iTotalGridSize = iGridSize1 * iGridSize2;
    dTimeStep = stod(argv[4]);
    dInitialEnergy = pow(10.0, stod(argv[5]));
    dDomainFactor1 = stod(argv[6]);
    dDomainFactor2 = stod(argv[7]);
    dDomainSize1 = 2 * dPI * stod(argv[6]);
    dDomainSize2 = 2 * dPI * stod(argv[7]);
    dTimeWindow = pow(10.0, stod(argv[8]));
    bOptimizeSolution = stoi(argv[9]);
    dOptimizationTolerance = stod(argv[10]);
    bNumericalContinuation = stoi(argv[11]);

    vLinearOperator.resize(iTotalGridSize);
    vLaplaceOperator.resize(iTotalGridSize);
    vDifferentialOperator1.resize(iTotalGridSize);
    vDifferentialOperator2.resize(iTotalGridSize);
    getMathematicalOperators();

    coeffAlphaI = { 343038331393.0 / 1130875731271.0,
                    288176579239.0 / 1140253497719.0,
                    253330171251.0 / 677500478386.0,
                    189462239225.0 / 1091147436423.0 };
    coeffBetaI = { 35965327958.0 / 140127563663.0,
                   19632212512.0 / 2700543775099.0,
                  -173747147147.0 / 351772688865.0,
                   91958533623.0 / 727726057489.0 };
    coeffAlphaE = { 14.0 / 25.0,
                    777974228744.0 / 1346157007247.0,
                    251277807242.0 / 1103637129625.0,
                    113091689455.0 / 220187950967.0 };
    coeffBetaE = { 0.0,
                  -251352885992.0 / 790610919619.0,
                  -383714262797.0 / 1103637129625.0,
                  -403360439203.0 / 1888264787188.0 };
    
    if (bOptimizeSolution == true)
        dOptimalTimeWindow = dTimeWindow;
    else
        dOptimalTimeWindow = pow(10.0, stod(argv[12]));

    cout << "Parameter settings:\nIC " << strInitialGuessName 
              << ", N_x1 " << iGridSize1 
              << ", N_x2 " << iGridSize2       
              << ", dt " << dTimeStep         
              << ", K " << dInitialEnergy    
              << ", ell1 " << dDomainFactor1     
              << ", ell2 " << dDomainFactor2  
              << ", L_1 " << dDomainSize1     
              << ", L_2 " << dDomainSize2     
              << ", T " << dTimeWindow       
              << ", nsteps " << iGetNumericalSteps() 
              << ", nstepsFile " << iGetNumericalStepsPerFile() 
              << ", opt " << bOptimizeSolution       
              << ", tol " << dOptimizationTolerance  
              << ", cont " << bNumericalContinuation  
              << ", optT " << dOptimalTimeWindow      
              << endl;
}

size_t Parameters::iGetNumericalSteps() const {   
    return static_cast<size_t>(round(dTimeWindow / dTimeStep));
}

size_t Parameters::iGetNumericalStepsPerFile() const {   
    return static_cast<size_t>(round( 4e7 / (iTotalGridSize) ));
}

void Parameters::getMathematicalOperators() {
    dEnergyFactor = dDomainFactor1 * dDomainFactor2 * (2.0 * dPI) * (2.0 * dPI) / (static_cast<double>(iTotalGridSize) * iTotalGridSize);
    
    _vGridpoints1 = setPhysicalSpace1();
    _vGridpoints2 = setPhysicalSpace2();

    const double _dSpaceStep1 = dDomainSize1 / static_cast<double>(iGridSize1);
    const double _dSpaceStep2 = dDomainSize2 / static_cast<double>(iGridSize2);
    dSpaceStep = min(_dSpaceStep1, _dSpaceStep2);

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

    vLinearOperator.resize(iTotalGridSize);
    vLaplaceOperator.resize(iTotalGridSize);
    vDifferentialOperator1.resize(iTotalGridSize);
    vDifferentialOperator2.resize(iTotalGridSize);
    vDealiasingOperator.resize(iTotalGridSize);
    const double kCut1 = (2.0 / 3.0) * (static_cast<double>(iGridSize1) / 2.0);
    const double kCut2 = (2.0 / 3.0) * (static_cast<double>(iGridSize2) / 2.0);

    double maxRadius = 0.0;
    vRadialBin.resize(iTotalGridSize);

    for (size_t i = 0; i < iGridSize2; ++i) {
        for (size_t j = 0; j < iGridSize1; ++j) {
            const size_t k = getIndex(i, j, iGridSize1);
            vGrid1[k] = _vGridpoints1[j];
            vGrid2[k] = _vGridpoints2[i];
            _vSpectrumNonlinear1[k] = vWavenumbersNonlinear1[j];
            _vSpectrumNonlinear2[k] = vWavenumbersNonlinear2[i];
            _vSpectrumLinear1[k] = _vWavenumbersLinear1[j];
            _vSpectrumLinear2[k] = _vWavenumbersLinear2[i];
            
            const double _dSpectralRadiusSquared = _vSpectrumLinear1[k] * _vSpectrumLinear1[k] + _vSpectrumLinear2[k] * _vSpectrumLinear2[k];
            const double _dBilaplacianValue = _dSpectralRadiusSquared * _dSpectralRadiusSquared;
            vH1Weight[k] = 1.0 + _dSpectralRadiusSquared;
            vH2Weight[k] = vH1Weight[k] * vH1Weight[k];
            
            vLinearOperator[k] = -_dSpectralRadiusSquared + _dBilaplacianValue;
            vLaplaceOperator[k] = -_dSpectralRadiusSquared;
            vDifferentialOperator1[k] = Imaginary * _vSpectrumNonlinear1[k];
            vDifferentialOperator2[k] = Imaginary * _vSpectrumNonlinear2[k];

            const bool keepMode = abs(vWavenumbersNonlinear1[j]) <= kCut1 && abs(vWavenumbersNonlinear2[i]) <= kCut2;
            vDealiasingOperator[k] = keepMode ? complex<double>{1.0, 0.0} : complex<double>{0.0, 0.0};

            const double k1 = _vWavenumbersLinear1[j];
            const double k2 = _vWavenumbersLinear2[i];
            const double radius = sqrt(k1 * k1 + k2 * k2);
            vRadialBin[k] = static_cast<size_t>(round(radius));
            maxRadius = max(maxRadius, radius);
        }   
    }
    iMaxRadialBin = static_cast<size_t>(ceil(maxRadius));
}