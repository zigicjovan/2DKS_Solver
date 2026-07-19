#include "Parameters.h"
#include <cmath> 
#include <string>

#include <iostream>
#include <vector>
#include <complex>
#include <algorithm>
#include <array>

using namespace std;

// Private functions
vector<double> Parameters::setPhysicalSpace1() {
    vector<double> gridpoints(_iGridSize1);

    for (size_t i = 0; i < _iGridSize1; ++i) {
        gridpoints[i] = _dDomainSize1 * static_cast<double>(i) / static_cast<double>(_iGridSize1);
    }

    return gridpoints;
}

vector<double> Parameters::setPhysicalSpace2() {
    vector<double> gridpoints(_iGridSize2);

    for (size_t i = 0; i < _iGridSize2; ++i) {
        gridpoints[i] = _dDomainSize2 * static_cast<double>(i) / static_cast<double>(_iGridSize2);
    }

    return gridpoints;
}

vector<double> Parameters::setFourierModesNonlinear1() {
    vector<double> wavenumbers(_iGridSize1);

    for (size_t i = 0; i < _iGridSize1 / 2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize1 * static_cast<double>(i);
    }

    wavenumbers[_iGridSize1 / 2] = 0.0;

    for (size_t i = _iGridSize1 / 2 + 1; i < _iGridSize1; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize1 * (static_cast<double>(i) - static_cast<double>(_iGridSize1));
    }

    return wavenumbers;
}

vector<double> Parameters::setFourierModesNonlinear2() {
    vector<double> wavenumbers(_iGridSize2);

    for (size_t i = 0; i < _iGridSize2 / 2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize2 * static_cast<double>(i);
    }

    wavenumbers[_iGridSize2 / 2] = 0.0;
    
    for (size_t i = _iGridSize2 / 2 + 1; i < _iGridSize2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize2 * (static_cast<double>(i) - static_cast<double>(_iGridSize2));
    }

    return wavenumbers;
}

vector<double> Parameters::setFourierModesLinear1() {
    vector<double> wavenumbers(_iGridSize1);
    
    for (size_t i = 0; i <= _iGridSize1 / 2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize1 * static_cast<double>(i);
    }

    for (size_t i = _iGridSize1 / 2 + 1; i < _iGridSize1; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize1 * (static_cast<double>(i) - static_cast<double>(_iGridSize1));
    }

    return wavenumbers;
}

vector<double> Parameters::setFourierModesLinear2() {
    vector<double> wavenumbers(_iGridSize2);
    
    for (size_t i = 0; i <= _iGridSize2 / 2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize2 * static_cast<double>(i);
    }
    
    for (size_t i = _iGridSize2 / 2 + 1; i < _iGridSize2; ++i) {
        wavenumbers[i] = 2.0 * dPI / _dDomainSize2 * (static_cast<double>(i) - static_cast<double>(_iGridSize2));
    }
    
    return wavenumbers;
}

// Public functions
Parameters::Parameters(int argc, char* argv[]) {

    _strInitialGuessName = argv[1];
    _iGridSize1 = stoi(argv[2]);
    _iGridSize2 = stoi(argv[3]);
    _iTotalGridSize = _iGridSize1 * _iGridSize2;
    _dTimeStep = stod(argv[4]);
    _dInitialEnergy = pow(10.0, stod(argv[5]));
    _dDomainFactor1 = stod(argv[6]);
    _dDomainFactor2 = stod(argv[7]);
    _dDomainSize1 = 2 * dPI * stod(argv[6]);
    _dDomainSize2 = 2 * dPI * stod(argv[7]);
    _dTimeWindow = pow(10.0, stod(argv[8]));
    _bOptimizeSolution = stoi(argv[9]);
    _dOptimizationTolerance = stod(argv[10]);
    _bNumericalContinuation = stoi(argv[11]);
    _bActiveLineSearch = false;

    if (_bOptimizeSolution) {
        _dOptimalTimeWindow = _dTimeWindow;
    }
    else {
        _dOptimalTimeWindow = pow(10.0, stod(argv[12]));
    }

    _iSavedStates = stoi(argv[13]);

    _vLinearOperator.resize(_iTotalGridSize);
    _vLaplaceOperator.resize(_iTotalGridSize);
    _vDifferentialOperator1.resize(_iTotalGridSize);
    _vDifferentialOperator2.resize(_iTotalGridSize);
    getMathematicalOperators();

    _coeffAlphaI = { 343038331393.0 / 1130875731271.0,
                    288176579239.0 / 1140253497719.0,
                    253330171251.0 / 677500478386.0,
                    189462239225.0 / 1091147436423.0 };
    _coeffBetaI = { 35965327958.0 / 140127563663.0,
                   19632212512.0 / 2700543775099.0,
                  -173747147147.0 / 351772688865.0,
                   91958533623.0 / 727726057489.0 };
    _coeffAlphaE = { 14.0 / 25.0,
                    777974228744.0 / 1346157007247.0,
                    251277807242.0 / 1103637129625.0,
                    113091689455.0 / 220187950967.0 };
    _coeffBetaE = { 0.0,
                  -251352885992.0 / 790610919619.0,
                  -383714262797.0 / 1103637129625.0,
                  -403360439203.0 / 1888264787188.0 };

    const int dRequiredMemory = 0.640 * (getNumericalSteps() + getNumericalStepsPerFile() - 1) / getNumericalStepsPerFile();
    const int dFinalMemory = _iSavedStates * dRequiredMemory / getNumericalSteps();

    cout << "Parameter settings:\nIC " << _strInitialGuessName 
              << ", N_x1 " << _iGridSize1 
              << ", N_x2 " << _iGridSize2       
              << ", dt " << _dTimeStep         
              << ", K " << _dInitialEnergy    
              << ", ell1 " << _dDomainFactor1     
              << ", ell2 " << _dDomainFactor2  
              << ", L_1 " << _dDomainSize1     
              << ", L_2 " << _dDomainSize2     
              << ", T " << _dTimeWindow       
              << ", opt " << _bOptimizeSolution       
              << ", tol " << _dOptimizationTolerance  
              << ", cont " << _bNumericalContinuation  
              << ", optT " << _dOptimalTimeWindow      
              << endl
              << "Intermediate Storage (for adjoint solve) = " << dRequiredMemory
              << " GB, Total Timesteps = " << getNumericalSteps() 
              << " (Max File Timesteps (640MB) " << getNumericalStepsPerFile() 
              << ", Remainder File Timesteps " << ( getNumericalSteps() % getNumericalStepsPerFile() )
              << ")" << endl
              << "Final Storage = " << dFinalMemory
              << " GB, Total Timesteps = " << _iSavedStates << endl;
}

void Parameters::getMathematicalOperators() {
    _dEnergyFactor = _dDomainFactor1 * _dDomainFactor2 * (2.0 * dPI) * (2.0 * dPI) / (static_cast<double>(_iTotalGridSize) * _iTotalGridSize);
    
    _vGridpoints1 = setPhysicalSpace1();
    _vGridpoints2 = setPhysicalSpace2();

    const double _dSpaceStep1 = _dDomainSize1 / static_cast<double>(_iGridSize1);
    const double _dSpaceStep2 = _dDomainSize2 / static_cast<double>(_iGridSize2);
    _dSpaceStep = min(_dSpaceStep1, _dSpaceStep2);

    _vWavenumbersNonlinear1 = setFourierModesNonlinear1();
    _vWavenumbersNonlinear2 = setFourierModesNonlinear2();

    _vWavenumbersLinear1 = setFourierModesLinear1();
    _vWavenumbersLinear2 = setFourierModesLinear2();

    _vGrid1.resize(_iTotalGridSize);
    _vGrid2.resize(_iTotalGridSize);

    _vSpectrumNonlinear1.resize(_iTotalGridSize);
    _vSpectrumNonlinear2.resize(_iTotalGridSize);
    _vSpectrumLinear1.resize(_iTotalGridSize);
    _vSpectrumLinear2.resize(_iTotalGridSize);

    _vH1Weight.resize(_iTotalGridSize);
    _vH2Weight.resize(_iTotalGridSize);

    _vLinearOperator.resize(_iTotalGridSize);
    _vLaplaceOperator.resize(_iTotalGridSize);
    _vDifferentialOperator1.resize(_iTotalGridSize);
    _vDifferentialOperator2.resize(_iTotalGridSize);
    _vDealiasingOperator.resize(_iTotalGridSize);
    const double kCut1 = (2.0 / 3.0) * (static_cast<double>(_iGridSize1) / 2.0);
    const double kCut2 = (2.0 / 3.0) * (static_cast<double>(_iGridSize2) / 2.0);

    double maxRadius = 0.0;
    _vRadialBin.resize(_iTotalGridSize);

    for (size_t i = 0; i < _iGridSize2; ++i) {
        for (size_t j = 0; j < _iGridSize1; ++j) {
            const size_t k = getIndex(i, j, _iGridSize1);
            _vGrid1[k] = _vGridpoints1[j];
            _vGrid2[k] = _vGridpoints2[i];
            _vSpectrumNonlinear1[k] = _vWavenumbersNonlinear1[j];
            _vSpectrumNonlinear2[k] = _vWavenumbersNonlinear2[i];
            _vSpectrumLinear1[k] = _vWavenumbersLinear1[j];
            _vSpectrumLinear2[k] = _vWavenumbersLinear2[i];
            
            const double _dSpectralRadiusSquared = _vSpectrumLinear1[k] * _vSpectrumLinear1[k] + _vSpectrumLinear2[k] * _vSpectrumLinear2[k];
            const double _dBilaplacianValue = _dSpectralRadiusSquared * _dSpectralRadiusSquared;
            _vH1Weight[k] = 1.0 + _dSpectralRadiusSquared;
            _vH2Weight[k] = _vH1Weight[k] * _vH1Weight[k];
            
            _vLinearOperator[k] = -_dSpectralRadiusSquared + _dBilaplacianValue;
            _vLaplaceOperator[k] = -_dSpectralRadiusSquared;
            _vDifferentialOperator1[k] = Imaginary * _vSpectrumNonlinear1[k];
            _vDifferentialOperator2[k] = Imaginary * _vSpectrumNonlinear2[k];

            const bool keepMode = abs(_vWavenumbersNonlinear1[j]) <= kCut1 && abs(_vWavenumbersNonlinear2[i]) <= kCut2;
            _vDealiasingOperator[k] = keepMode ? complex<double>{1.0, 0.0} : complex<double>{0.0, 0.0};

            const double k1 = _vWavenumbersLinear1[j];
            const double k2 = _vWavenumbersLinear2[i];
            const double radius = sqrt(k1 * k1 + k2 * k2);
            _vRadialBin[k] = static_cast<size_t>(round(radius));
            maxRadius = max(maxRadius, radius);
        }   
    }
    _iMaxRadialBin = static_cast<size_t>(ceil(maxRadius));
}

size_t Parameters::getIndex(size_t i, size_t j, size_t N) {
    return i * N + j;
}

size_t Parameters::getNumericalSteps() const {   
    return static_cast<size_t>(round(_dTimeWindow / _dTimeStep));
}

size_t Parameters::getNumericalStepsPerFile() const {   
    return static_cast<size_t>(round( 4e7 / (_iTotalGridSize) ));
}

const vector<double>& Parameters::getGrid1() const { 
    return _vGrid1; 
}

const vector<double>& Parameters::getGrid2() const { 
    return _vGrid2; 
}

const vector<double>& Parameters::getWavenumbersNonlinear1() const { 
    return _vWavenumbersNonlinear1; 
}

const vector<double>& Parameters::getWavenumbersNonlinear2() const { 
    return _vWavenumbersNonlinear2; 
}

const string& Parameters::getInitialGuessName() const { 
    return _strInitialGuessName; 
}

size_t Parameters::getGridSize1() const { 
    return _iGridSize1; 
}

size_t Parameters::getGridSize2() const { 
    return _iGridSize2; 
}

size_t Parameters::getTotalGridSize() const { 
    return _iTotalGridSize; 
}

size_t Parameters::getSavedStates() const { 
    return _iSavedStates; 
}

double Parameters::getTimeStep() const { 
    return _dTimeStep; 
}

double Parameters::getSpaceStep() const { 
    return _dSpaceStep; 
}

double Parameters::getInitialEnergy() const { 
    return _dInitialEnergy; 
}

double Parameters::getEnergyFactor() const { 
    return _dEnergyFactor; 
}

double Parameters::getDomainFactor1() const { 
    return _dDomainFactor1; 
}

double Parameters::getDomainFactor2() const { 
    return _dDomainFactor2; 
}

double Parameters::getDomainSize1() const { 
    return _dDomainSize1; 
}

double Parameters::getDomainSize2() const { 
    return _dDomainSize2; 
}

double Parameters::getTimeWindow() const { 
    return _dTimeWindow; 
}

const vector<complex<double>>& Parameters::getLinearOperator() const {
    return _vLinearOperator;
}

const vector<complex<double>>& Parameters::getLaplaceOperator() const {
    return _vLaplaceOperator;
}

const vector<complex<double>>& Parameters::getDifferentialOperator1() const {
    return _vDifferentialOperator1;
}

const vector<complex<double>>& Parameters::getDifferentialOperator2() const {
    return _vDifferentialOperator2;
}

const vector<complex<double>>& Parameters::getDealiasingOperator() const {
    return _vDealiasingOperator;
}

const vector<double>& Parameters::getH1Weight() const { 
    return _vH1Weight; 
}

const vector<double>& Parameters::getH2Weight() const { 
    return _vH2Weight; 
}

const vector<size_t>& Parameters::getRadialBin() const { 
    return _vRadialBin; 
}

size_t Parameters::getMaxRadialBin() const { 
    return _iMaxRadialBin; 
}

const array<double, 4>& Parameters::getCoeffAlphaI() const { 
    return _coeffAlphaI; 
}

const array<double, 4>& Parameters::getCoeffBetaI() const { 
    return _coeffBetaI; 
}

const array<double, 4>& Parameters::getCoeffAlphaE() const { 
    return _coeffAlphaE; 
}

const array<double, 4>& Parameters::getCoeffBetaE() const { 
    return _coeffBetaE; 
}

bool Parameters::getOptimizeSolution() const { 
    return _bOptimizeSolution; 
}

void Parameters::setOptimizeSolution(bool optimizeSolution) {
    _bOptimizeSolution = optimizeSolution;
}

bool Parameters::getActiveLineSearch() const { 
    return _bActiveLineSearch; 
}

void Parameters::setActiveLineSearch(bool lineSearch) {
    _bActiveLineSearch = lineSearch;
}

double Parameters::getOptimizationTolerance() const {
    return _dOptimizationTolerance;
}

bool Parameters::getNumericalContinuation() const {
    return _bNumericalContinuation;
}

double Parameters::getOptimalTimeWindow() const {
    return _dOptimalTimeWindow;
}