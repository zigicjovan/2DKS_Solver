#include "Solver.h"
#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "FFTWPlanner.h"
#include "Timer.h"

#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <array>
#include <vector>

#include <filesystem>
#include <string>
#include <iostream>

#include <complex>
#include <cstddef>
#include <random>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;

// TODO: create standalone functions for big cases

// Private functions
void Solver::saveSolutionDiagnostics(const vector<array<double, 4>>& vDiagnostics) {
    ofstream file(_paths.fEnergyEvolution);
    file << setprecision(16) << scientific;
    for (const auto& row : vDiagnostics)
        file << row[0] << ' ' << row[1] << ' ' << row[2] << ' ' << row[3] << '\n';
}

void Solver::saveSolutionSpectrum(const vector<vector<double>>& vSpectrumHistory) {
    ofstream file(_paths.fFourierSpectrumEvolution);
    file << setprecision(16) << scientific;
    for (size_t i = 0; i <= _params.iMaxRadialBin; ++i) {
        file << i;
        for (const auto& spectrum : vSpectrumHistory)
            file << ' ' << spectrum[i];
        file << '\n';
    }
}

void Solver::checkCFL(SolutionData& vData1, SolutionData& vData2) {
    constexpr double CFL_LIMIT = 0.9;
    double dMaxGradientMagnitude = 0.0;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i) {
        const double vel = hypot(vData1[i].real(), vData2[i].real()); // sqrt(ux * ux + uy * uy)
        dMaxGradientMagnitude = max(dMaxGradientMagnitude, vel);
    }
    const double cflValue = _params.dTimeStep * dMaxGradientMagnitude / _params.dSpaceStep;
    if (cflValue > CFL_LIMIT) {
        cerr << "ERROR: CFL condition violated." << " Maximum gradient magnitude = " << dMaxGradientMagnitude << ", CFL value = " << cflValue << ", CFL Limit = " << CFL_LIMIT << '\n';
        exit(EXIT_FAILURE);
    }
}

// Public functions
Solver::Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer) 
                : _params(params), _paths(paths), _fftwPlan(fftwPlan), _timer(timer) {};

void Solver::setInitialCondition(SolutionData& vTargetState) {  
    if (_params.strInitialGuessName == "s1") {
        size_t stateNumber = 0;
        for (size_t i = 0; i < _params.iGridSize2; ++i) {
            for (size_t j = 0; j < _params.iGridSize1; ++j) {
                const size_t k = getIndex(i, j, _params.iGridSize1);
                const double x = _params.vGrid1[k];
                const double y = _params.vGrid2[k];
                vTargetState(j, i, stateNumber) = complex<double>(
                    sin(x / _params.dDomainFactor1 + y / _params.dDomainFactor2)
                    + sin(x / _params.dDomainFactor1) + sin(y / _params.dDomainFactor2), 0.0);
            }
        }
        _fftwPlan.fft2InPlace(vTargetState.getDataPointer());
    }
    else if (_params.strInitialGuessName == "randfour") {
        const size_t stateNumber = 0;

        const double kmax  = max(_params.iGridSize1, _params.iGridSize2) / 2.0;
        const double delta = 6.0 / kmax;

        mt19937 rng(random_device{}());
        uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);

        complex<double>* vTempState = vTargetState.getStatePointer(stateNumber);

        // Build randomized Fourier coefficients
        for (size_t i = 0; i < _params.iGridSize2; ++i) {
            for (size_t j = 0; j < _params.iGridSize1; ++j) {
                const size_t k = getIndex(i, j, _params.iGridSize1);
                const double k1 = _params.vWavenumbersNonlinear1[j];
                const double k2 = _params.vWavenumbersNonlinear2[i];
                const double kabs = sqrt(k1 * k1 + k2 * k2);
                const double A = exp(-delta * kabs);
                const double phi = dist(rng);
                vTempState[k] = A * complex<double>(cos(phi), sin(phi));
            }
        }
        vTempState[0] = complex<double>(0.0, 0.0); // zero mean mode
        _fftwPlan.ifft2InPlace(_params, vTempState);  // real(ifft2(uhat))
        for (size_t i = 0; i < _params.iTotalGridSize; ++i)
            vTempState[i] = complex<double>(vTempState[i].real(), 0.0);
        _fftwPlan.fft2InPlace(vTempState);
    }
    else {
        throw runtime_error(
            "Unknown initial condition: " + _params.strInitialGuessName);
    }
    vTargetState.setInitialEnergyL2();
}

void Solver::setSolutionState(StateSolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        
        case SolveInitialState: {               
            if (_params.bNumericalContinuation == 0) {
                setInitialCondition(vTargetState);
                cout << "Generated initial guess.\n";
            }
            else if (filesystem::exists(_paths.fOptimalInitialData)) {
                vTargetState.loadData(InitialState);
                cout << "Loaded optimal initial data.\n";
            }
            else {
                filesystem::path fContinuedFile;
                filesystem::path fTempFile;
                double dBestT = _params.dTimeWindow + 1.0;            
                
                // Search for nearest previous time window in directory
                for (const auto& entry : filesystem::directory_iterator(_paths.dirOptimalInitialData)) {
                    if (!entry.is_regular_file())
                        continue;
                    string filename = entry.path().filename().string();
                    if (filename.find(_paths.strTestcaseGeneric.str()) == string::npos)
                        continue;
                    size_t pos1 = filename.find("_T_"); 
                    if (pos1 == string::npos)
                        continue;
                    pos1 += 3; // skip "_T_"
                    size_t pos2 = filename.find("_opt_", pos1); 
                    if (pos2 == string::npos)
                        continue;
                    double T = stod(filename.substr(pos1, pos2 - pos1)); // Extract T
                    if (T <= _params.dTimeWindow && ( T > dBestT || dBestT > _params.dTimeWindow)) {
                        dBestT = T;
                        fContinuedFile = entry.path();
                    }
                }
                if (dBestT <= _params.dTimeWindow) {
                    fTempFile = _paths.fOptimalInitialData;
                    _paths.fOptimalInitialData = fContinuedFile;
                    vTargetState.loadData(InitialState);
                    cout << "Loaded continued initial data from T = " << dBestT << "\n" ;
                    _paths.fOptimalInitialData = fTempFile; 
                }
                else {
                    setInitialCondition(vTargetState);
                    cout << "No saved data, generated initial guess.\n";
                }
            }
            vTargetState.saveData(InitialState);
            break;
        }

        case SolveTerminalState: {
            vTargetState.saveData(TerminalState);
            break;
        }
    }  
}

void Solver::setSolutionInTime(TimeSolutionType targetType, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        
        case SolveForwardInTime: {
            double dTimePoint = 0.0;
            const size_t totalSteps = _params.iGetNumericalSteps();
            const size_t stepsPerFile = _params.iGetNumericalStepsPerFile();
            const size_t remainderSteps = totalSteps % stepsPerFile;
            const size_t fullSteps = totalSteps - remainderSteps;

            vector<array<double, 4>> vDiagnostics;
            vDiagnostics.reserve(_params.iGetNumericalSteps() + 1);
            vector<vector<double>> vSpectrumHistory;
            vSpectrumHistory.reserve(_params.iGetNumericalSteps() + 1);

            SolutionData vStateCurrent = vTargetStart; // phi_{i} = phi_{0}
            SolutionData vStateNext(_params, _paths, InitialState); // phi_{i+1}
            SolutionData vNonlinearTermPrevious(_params, _paths, InitialState); // N(phi_{i-1})
            for (size_t p = 0; p < _params.iTotalGridSize; ++p)
                vNonlinearTermPrevious[p] = complex<double>{0.0, 0.0};
            SolutionData vNonlinearTermCurrent(_params, _paths, InitialState); // N(phi_{i})
            SolutionData vFunctionSpatialDerivative1(_params, _paths, InitialState); // (phi_{i})_{x_1}
            SolutionData vFunctionSpatialDerivative2(_params, _paths, InitialState); // (phi_{i})_{x_2}

            vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2() });
            vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
            for (size_t i = 1; i < _params.iGetNumericalSteps() + 1; ++i) {
                dTimePoint = i * _params.dTimeStep;
                for (size_t j = 0; j < 4; ++j) {
                    
                    // f_x and f_y in Fourier space
                    for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                        vFunctionSpatialDerivative1[k] = _params.vDifferentialOperator1[k] * vStateCurrent[k];
                        vFunctionSpatialDerivative2[k] = _params.vDifferentialOperator2[k] * vStateCurrent[k];
                    }
                    
                    // pseudospectral products f_x^2 and f_y^2 in physical space
                    _fftwPlan.ifft2InPlace(_params, vFunctionSpatialDerivative1.getDataPointer());
                    _fftwPlan.ifft2InPlace(_params, vFunctionSpatialDerivative2.getDataPointer());
                    checkCFL(vFunctionSpatialDerivative1, vFunctionSpatialDerivative2);

                    // build nonlinear term in physical space
                    for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                        const double ux = vFunctionSpatialDerivative1[k].real();
                        const double uy = vFunctionSpatialDerivative2[k].real();
                        vNonlinearTermCurrent[k] = complex<double>{0.5 * (ux * ux + uy * uy), 0.0};
                    }

                    // dealias transformed nonlinear term and IMEX RK4 update
                    _fftwPlan.fft2InPlace(vNonlinearTermCurrent.getDataPointer());
                    for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                        const complex<double> Lin = _params.vLinearOperator[k];
                        vNonlinearTermCurrent[k] *= _params.vDealiasingOperator[k];
                        vStateNext[k] = ( (1.0 - _params.dTimeStep * _params.coeffBetaI[j] * Lin) * vStateCurrent[k]
                            - _params.dTimeStep * _params.coeffAlphaE[j] * vNonlinearTermCurrent[k]
                            - _params.dTimeStep * _params.coeffBetaE[j]  * vNonlinearTermPrevious[k] ) 
                            / (1.0 + _params.dTimeStep * _params.coeffAlphaI[j] * Lin);
                    }
                    vStateCurrent.swapDataFrom(vStateNext);
                    vNonlinearTermPrevious.swapDataFrom(vNonlinearTermCurrent);
                }
                vStateCurrent[0] = complex<double>{0.0, 0.0}; // mean-zero correction

                if (i <= fullSteps) {
                    const size_t localIndex = (i - 1) % stepsPerFile;
                    vHistoryIntermediate.setData(vStateCurrent, localIndex);
                    if (i % stepsPerFile == 0)
                        vHistoryIntermediate.saveData(IntermediateHistory, dTimePoint);
                }
                else {
                    const size_t localIndex = i - fullSteps - 1;
                    vHistoryRemainder.setData(vStateCurrent, localIndex);
                    if (i == totalSteps)
                        vHistoryRemainder.saveData(RemainderHistory, dTimePoint);
                }

                vDiagnostics.push_back({dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2()});
                vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
            }
            vTargetEnd.moveDataFrom(vStateCurrent);
            setSolutionState(SolveTerminalState, vTargetEnd);
            _timer.printInterval("Forward problem solved at ");
            if (_params.bOptimizeSolution == 0) {
                saveSolutionDiagnostics(vDiagnostics); 
                saveSolutionSpectrum(vSpectrumHistory);
            }
            cout << "\n";
            break;
        }

        case SolveBackwardInTime: {
            
            // TO DO: solve backward problem using IMEX RK4 [create function]
            // loadData(_paths, vHistoryIntermediate, IntermediateHistory, _params.dTimeWindow);
            // loadData(_paths, vHistoryRemainder, RemainderHistory, dCurrentT);
            // saveData: BackwardInitialState
            _timer.printInterval("Backward problem solved at ");
            cout << "\n";
            break;
        }
    }  
}

double Solver::getOptimalSolution(OptimizeSolutionType targetType, SolutionData& vObjectiveGradient, SolutionData& vTargetStart, 
                SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    double targetValue = 0.0; // energy or step size
    switch (targetType) {
        
        case OptimizeEnergyAmplification: {
            if (_params.bOptimizeSolution == 1) {
                
                // TO DO: solve RCG  [create function]
                // saveData(vTargetStart, InitialState);
                // save SolutionBranch, OptDiagnostics
                
                // save diagnostics after completion
                _params.bOptimizeSolution = 0; 
                setSolutionInTime(SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            }
            _timer.printInterval("Energy maximization problem solved at ");
            cout << "\n";
            break;
        }

        case OptimizeLineSearchStepSize: {
            
            // TO DO: solve Brent [create function]
            // save OptLineSearch
            _timer.printInterval("Step size optimization problem solved at ");
            cout << "\n";
            break;
        }
    } 
    return targetValue;
}
