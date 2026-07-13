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

// Index for access, push_back for assignment
// TODO: create standalone functions for big cases

// Private functions
void Solver::saveSolutionDiagnostics(const vector<array<double, 4>>& vDiagnostics) {
    ofstream file(paths.fEnergyEvolution);
    file << setprecision(16) << scientific;
    for (const auto& row : vDiagnostics)
        file << row[0] << ' ' << row[1] << ' ' << row[2] << ' ' << row[3] << '\n';
}

void Solver::saveSolutionSpectrum(const vector<vector<double>>& vSpectrumHistory) {
    ofstream file(paths.fFourierSpectrumEvolution);
    file << setprecision(16) << scientific;
    for (size_t i = 0; i <= params.iMaxRadialBin; ++i) {
        file << i;
        for (const auto& spectrum : vSpectrumHistory)
            file << ' ' << spectrum[i];
        file << '\n';
    }
}

void Solver::checkCFL(SolutionData& vData1, SolutionData& vData2) {
    constexpr double CFL_LIMIT = 0.9;
    double dMaxGradientMagnitude = 0.0;
    for (size_t i = 0; i < params.iTotalGridSize; ++i) {
        const double vel = hypot(vData1[i].real(), vData2[i].real()); // sqrt(ux * ux + uy * uy)
        dMaxGradientMagnitude = max(dMaxGradientMagnitude, vel);
    }
    const double cflValue = params.dTimeStep * dMaxGradientMagnitude / params.dSpaceStep;
    if (cflValue > CFL_LIMIT) {
        cerr << "ERROR: CFL condition violated." << " Maximum gradient magnitude = " << dMaxGradientMagnitude << ", CFL value = " << cflValue << ", CFL Limit = " << CFL_LIMIT << '\n';
        exit(EXIT_FAILURE);
    }
}

// Public functions
Solver::Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer) 
                : params(params), paths(paths), fftwPlan(fftwPlan), timer(timer) {};

void Solver::setInitialCondition(SolutionData& vTargetState) {  
    if (params.strInitialGuessName == "s1") {
        int stateNumber = 0;
        for (size_t i = 0; i < params.iGridSize2; ++i) {
            for (size_t j = 0; j < params.iGridSize1; ++j) {
                const size_t k = getIndex(i, j, params.iGridSize1);
                const double x = params.vGrid1[k];
                const double y = params.vGrid2[k];
                vTargetState(params, j, i, stateNumber) = complex<double>(
                    sin(x / params.dDomainFactor1 + y / params.dDomainFactor2)
                    + sin(x / params.dDomainFactor1) + sin(y / params.dDomainFactor2), 0.0);
            }
        }
        fftwPlan.fft2InPlace(vTargetState.getDataPointer());
    }
    else if (params.strInitialGuessName == "randfour") {
        const int stateNumber = 0;

        const double kmax  = max(params.iGridSize1, params.iGridSize2) / 2.0;
        const double delta = 6.0 / kmax;

        mt19937 rng(random_device{}());
        uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);

        complex<double>* vTempState = vTargetState.getStatePointer(params, stateNumber);

        // Build randomized Fourier coefficients
        for (size_t i = 0; i < params.iGridSize2; ++i) {
            for (size_t j = 0; j < params.iGridSize1; ++j) {
                const size_t k = getIndex(i, j, params.iGridSize1);
                const double k1 = params.vWavenumbersNonlinear1[j];
                const double k2 = params.vWavenumbersNonlinear2[i];
                const double kabs = sqrt(k1 * k1 + k2 * k2);
                const double A = exp(-delta * kabs);
                const double phi = dist(rng);
                vTempState[k] = A * complex<double>(cos(phi), sin(phi));
            }
        }
        vTempState[0] = complex<double>(0.0, 0.0); // zero mean mode
        fftwPlan.ifft2InPlace(params, vTempState);  // real(ifft2(uhat))
        for (size_t i = 0; i < params.iTotalGridSize; ++i)
            vTempState[i] = complex<double>(vTempState[i].real(), 0.0);
        fftwPlan.fft2InPlace(vTempState);
    }
    else {
        throw runtime_error(
            "Unknown initial condition: " + params.strInitialGuessName);
    }
    vTargetState.setInitialEnergyL2(params);
}

void Solver::setSolutionState(StateSolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        
        case SolveInitialState: {               
            if (params.bNumericalContinuation == 0) {
                setInitialCondition(vTargetState);
                cout << "Generated initial guess.\n";
            }
            else if (filesystem::exists(paths.fOptimalInitialData)) {
                vTargetState.loadData(paths, InitialState);
                cout << "Loaded optimal initial data.\n";
            }
            else {
                filesystem::path fContinuedFile;
                filesystem::path fTempFile;
                double dBestT = params.dTimeWindow + 1.0;            
                
                // Search for nearest previous time window in directory
                for (const auto& entry : filesystem::directory_iterator(paths.dirOptimalInitialData)) {
                    if (!entry.is_regular_file())
                        continue;
                    string filename = entry.path().filename().string();
                    if (filename.find(paths.strTestcaseGeneric.str()) == string::npos)
                        continue;
                    size_t pos1 = filename.find("_T_"); 
                    if (pos1 == string::npos)
                        continue;
                    pos1 += 3; // skip "_T_"
                    size_t pos2 = filename.find("_opt_", pos1); 
                    if (pos2 == string::npos)
                        continue;
                    double T = stod(filename.substr(pos1, pos2 - pos1)); // Extract T
                    if (T <= params.dTimeWindow && ( T > dBestT || dBestT > params.dTimeWindow)) {
                        dBestT = T;
                        fContinuedFile = entry.path();
                    }
                }
                if (dBestT <= params.dTimeWindow) {
                    fTempFile = paths.fOptimalInitialData;
                    paths.fOptimalInitialData = fContinuedFile;
                    vTargetState.loadData(paths, InitialState);
                    cout << "Loaded continued initial data from T = " << dBestT << "\n" ;
                    paths.fOptimalInitialData = fTempFile; 
                }
                else {
                    setInitialCondition(vTargetState);
                    cout << "No saved data, generated initial guess.\n";
                }
            }
            vTargetState.saveData(paths, InitialState);
            break;
        }

        case SolveTerminalState: {
            
            // TO DO: take final state from vHistoryRemainder
            vTargetState.saveData(paths, TerminalState);
            break;
        }
    }  
}

void Solver::setSolutionInTime(TimeSolutionType targetType, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        
        case SolveForwardInTime: {
            double dTimePoint = 0.0;
            vector<array<double, 4>> vDiagnostics;
            vDiagnostics.reserve(params.iGetNumericalSteps());
            vector<vector<double>> vSpectrumHistory;
            vSpectrumHistory.reserve(params.iGetNumericalSteps());

            SolutionData vStateCurrent = vTargetStart; // phi_{i} = phi_{0}
            SolutionData vStateNext(params, InitialState); // phi_{i+1}
            SolutionData vNonlinearTermPrevious(params, InitialState); // N(phi_{i-1})
            for (size_t p = 0; p < params.iTotalGridSize; ++p)
                vNonlinearTermPrevious[p] = complex<double>{0.0, 0.0};
            SolutionData vNonlinearTermCurrent(params, InitialState); // N(phi_{i})
            SolutionData vFunctionSpatialDerivative1(params, InitialState); // (phi_{i})_{x_1}
            SolutionData vFunctionSpatialDerivative2(params, InitialState); // (phi_{i})_{x_2}

            vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(params), vStateCurrent.getEnergyH1(params), vStateCurrent.getEnergyH2(params) });
            vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum(params));
            for (size_t i = 1; i < params.iGetNumericalSteps(); ++i) {
                dTimePoint = i * params.dTimeStep;
                for (size_t j = 0; j < 4; ++j) {
                    
                    // f_x and f_y in Fourier space
                    for (size_t k = 0; k < params.iTotalGridSize; ++k) {
                        vFunctionSpatialDerivative1[k] = params.vDifferentialOperator1[k] * vStateCurrent[k];
                        vFunctionSpatialDerivative2[k] = params.vDifferentialOperator2[k] * vStateCurrent[k];
                    }
                    
                    // pseudospectral products f_x^2 and f_y^2 in physical space
                    fftwPlan.ifft2InPlace(params, vFunctionSpatialDerivative1.getDataPointer());
                    fftwPlan.ifft2InPlace(params, vFunctionSpatialDerivative2.getDataPointer());
                    checkCFL(vFunctionSpatialDerivative1, vFunctionSpatialDerivative2);

                    // build nonlinear term in physical space
                    for (size_t k = 0; k < params.iTotalGridSize; ++k) {
                        const double ux = vFunctionSpatialDerivative1[k].real();
                        const double uy = vFunctionSpatialDerivative2[k].real();
                        vNonlinearTermCurrent[k] = complex<double>{0.5 * (ux * ux + uy * uy), 0.0};
                    }

                    // dealias transformed nonlinear term and IMEX RK4 update
                    fftwPlan.fft2InPlace(vNonlinearTermCurrent.getDataPointer());
                    for (size_t k = 0; k < params.iTotalGridSize; ++k) {
                        const complex<double> Lin = params.vLinearOperator[k];
                        vNonlinearTermCurrent[k] *= params.vDealiasingOperator[k];
                        vStateNext[k] = ( (1.0 - params.dTimeStep * params.coeffBetaI[j] * Lin) * vStateCurrent[k]
                            - params.dTimeStep * params.coeffAlphaE[j] * vNonlinearTermCurrent[k]
                            - params.dTimeStep * params.coeffBetaE[j]  * vNonlinearTermPrevious[k] ) 
                            / (1.0 + params.dTimeStep * params.coeffAlphaI[j] * Lin);
                    }
                    swap(vStateCurrent, vStateNext);
                    swap(vNonlinearTermPrevious, vNonlinearTermCurrent);
                }
                vStateCurrent[0] = complex<double>{0.0, 0.0}; // mean-zero correction
                // saveData: IntermediateHistory, RemainderHistory, TerminalState
                // saveData(paths, vHistoryIntermediate, IntermediateHistory, params.dTimeWindow);
                // saveData(paths, vHistoryRemainder, RemainderHistory, dCurrentT);
                vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(params), vStateCurrent.getEnergyH1(params), vStateCurrent.getEnergyH2(params) });
                vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum(params));
            }
            timer.printInterval("Forward problem solved at ");
            if (params.bOptimizeSolution == 0) {
                saveSolutionDiagnostics(vDiagnostics); 
                saveSolutionSpectrum(vSpectrumHistory);
            }
            cout << "\n";
            break;
        }

        case SolveBackwardInTime: {
            
            // TO DO: solve backward problem using IMEX RK4 [create function]
            // loadData(paths, vHistoryIntermediate, IntermediateHistory, params.dTimeWindow);
            // loadData(paths, vHistoryRemainder, RemainderHistory, dCurrentT);
            // saveData: BackwardInitialState
            timer.printInterval("Backward problem solved at ");
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
            if (params.bOptimizeSolution == 1) {
                
                // TO DO: solve RCG  [create function]
                // saveData(paths, vTargetStart, InitialState);
                // save SolutionBranch, OptDiagnostics
                
                // save diagnostics after completion
                params.bOptimizeSolution = 0; 
                setSolutionInTime(SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            }
            timer.printInterval("Energy maximization problem solved at ");
            cout << "\n";
            break;
        }

        case OptimizeLineSearchStepSize: {
            
            // TO DO: solve Brent [create function]
            // save OptLineSearch
            timer.printInterval("Step size optimization problem solved at ");
            cout << "\n";
            break;
        }
    } 
    return targetValue;
}
