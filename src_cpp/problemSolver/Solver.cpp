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
#include <limits>
#include <utility>

using namespace std;

// TO DO LATER: Add MPI

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

void Solver::saveOptimizationDiagnostics(const array<double, 7>& vDiagnostics) {
    ofstream file(_paths.fOptimizationDiagnostics, ios::app);
    file << setprecision(16) << scientific;
    file << vDiagnostics[0] << ' ' << vDiagnostics[1] << ' ' << vDiagnostics[2] << ' ' << vDiagnostics[3] 
         << ' ' << vDiagnostics[4] << ' ' << vDiagnostics[5] << ' ' << vDiagnostics[6] << '\n';
}

void Solver::saveSolutionBranchAndPowerLaws(const array<double, 7>& vOptimalEnergySolution) {
    vector<array<double, 7>> vBranch;
    vector<array<double, 7>> vK;
    vector<array<double, 7>> vL;
    vector<array<double, 7>> vTK;
    vector<array<double, 7>> vTL;

    {
        ifstream inputBranch(_paths.fOptimalSolutionBranch);
        ifstream inputK(_paths.fOptimalSolutionInitialEnergyPowerLaw);
        ifstream inputL(_paths.fOptimalSolutionDomainSizePowerLaw);
        ifstream inputTK(_paths.fOptimalSolutionEnergyTimeWindowPowerLaw);
        ifstream inputTL(_paths.fOptimalSolutionDomainTimeWindowPowerLaw);
        array<double, 7> solutionBranch;
        array<double, 7> solutionK;
        array<double, 7> solutionL;
        array<double, 7> solutionTK;
        array<double, 7> solutionTL;
        while (inputBranch >> solutionBranch[0] >> solutionBranch[1] >> solutionBranch[2] >> solutionBranch[3] >> solutionBranch[4] >> solutionBranch[5] >> solutionBranch[6])
            vBranch.push_back(solutionBranch);
        while (inputK >> solutionK[0] >> solutionK[1] >> solutionK[2] >> solutionK[3] >> solutionK[4] >> solutionK[5] >> solutionK[6])
            vK.push_back(solutionK);
        while (inputL >> solutionL[0] >> solutionL[1] >> solutionL[2] >> solutionL[3] >> solutionL[4] >> solutionL[5] >> solutionL[6])
            vL.push_back(solutionL);
        while (inputTK >> solutionTK[0] >> solutionTK[1] >> solutionTK[2] >> solutionTK[3] >> solutionTK[4] >> solutionTK[5] >> solutionTK[6])
            vTK.push_back(solutionTK);
        while (inputTL >> solutionTL[0] >> solutionTL[1] >> solutionTL[2] >> solutionTL[3] >> solutionTL[4] >> solutionTL[5] >> solutionTL[6])
            vTL.push_back(solutionTL);
    }

    auto updateAndSave = [&](vector<array<double, 7>>& solutions, const string& filename, auto matches, auto sortingRule) {
        auto existing = find_if(solutions.begin(), solutions.end(), [&](const auto& solution) { return matches(solution, vOptimalEnergySolution);});
        if (existing == solutions.end())
            solutions.push_back(vOptimalEnergySolution);
        else if (vOptimalEnergySolution[6] > (*existing)[6])
            *existing = vOptimalEnergySolution;

        sort(solutions.begin(), solutions.end(), sortingRule);
        ofstream output(filename);
        output << scientific << setprecision(16);
        for (const auto& solution : solutions) {
            for (size_t column = 0; column < solution.size(); ++column) {
                if (column > 0)
                    output << ' ';
                output << solution[column];
            }
            output << '\n';
        }
    };

    updateAndSave(vBranch, _paths.fOptimalSolutionBranch,
        [](const auto& a, const auto& b) {
            return a[3] == b[3];
        },
        [](const auto& a, const auto& b) {
            return a[3] < b[3];
        });

    updateAndSave(vK, _paths.fOptimalSolutionInitialEnergyPowerLaw,
        [](const auto& a, const auto& b) {
            return a[0] == b[0];
        },
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    updateAndSave(vL, _paths.fOptimalSolutionDomainSizePowerLaw,
        [](const auto& a, const auto& b) {
            return a[1] == b[1] && a[2] == b[2];
        },
        [](const auto& a, const auto& b) {
            if (a[1] != b[1])
                return a[1] < b[1];
            return a[2] < b[2];
        });

    updateAndSave(vTK, _paths.fOptimalSolutionEnergyTimeWindowPowerLaw,
        [](const auto& a, const auto& b) {
            return a[0] == b[0];
        },
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    updateAndSave(vTL, _paths.fOptimalSolutionDomainTimeWindowPowerLaw,
        [](const auto& a, const auto& b) {
            return a[1] == b[1] && a[2] == b[2];
        },
        [](const auto& a, const auto& b) {
            if (a[1] != b[1])
                return a[1] < b[1];
            return a[2] < b[2];
        });
}

void Solver::saveLineSearch(vector<double>& vLineSearchHistory) {
    ofstream file(_paths.fOptimizationLineSearch, ios::app);
    const double dNaN = numeric_limits<double>::quiet_NaN();
    file << setprecision(16) << scientific;
    for (size_t i = 0; i < 1000; ++i) {
        if (i < vLineSearchHistory.size())
            file << vLineSearchHistory[i];
        else
            file << dNaN;

        if (i + 1 < 1000)
            file << ' ';
    }
    file << '\n';
}

EnergyData Solver::getMaxEnergyL2InTimeWindow() {
    ifstream inputFile(_paths.fEnergyEvolution);
    double dTimepoint;
    double dEnergyL2;
    double dEnergyH1;
    double dEnergyH2;
    EnergyData result{ -numeric_limits<double>::infinity(), 0.0 };
    while (inputFile >> dTimepoint >> dEnergyL2 >> dEnergyH1 >> dEnergyH2) {
        if (dEnergyL2 > result.dEnergy) {
            result.dEnergy = dEnergyL2;
            result.dTimepoint = dTimepoint;
        }
    }
    return result;
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

void Solver::setInitialCondition(SolutionData& vTargetState) {  
    if (_params.strInitialGuessName == "s1") {
        size_t stateNumber = 0;
        for (size_t i = 0; i < _params.iGridSize2; ++i) {
            for (size_t j = 0; j < _params.iGridSize1; ++j) {
                const size_t k = _params.getIndex(i, j, _params.iGridSize1);
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
                const size_t k = _params.getIndex(i, j, _params.iGridSize1);
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
        cerr << "ERROR: Unknown initial condition: " << _params.strInitialGuessName << endl;
        exit(EXIT_FAILURE);
    }
    vTargetState.setInitialEnergyL2();
}

void Solver::findContinuationForInitialData(SolutionData& vTargetState) {
    filesystem::path fContinuedFile;
    filesystem::path fTempFile;
    double dBestT = _params.dTimeWindow + 1.0;            

    // Search for nearest previous time window in directory
    for (const auto& entry : filesystem::directory_iterator(_paths.dirOptimalInitialData)) {
        
        if (!entry.is_regular_file())
            continue;
        
        string filename = entry.path().filename().string();
        if (filename.find(_paths.strTestcaseGenericTime.str()) == string::npos)
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

void Solver::solveForwardInTime(SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    double dTimePoint = 0.0;
    const size_t totalSteps = _params.iGetNumericalSteps();
    const size_t stepsPerFile = _params.iGetNumericalStepsPerFile();
    const size_t remainderSteps = totalSteps % stepsPerFile;
    const size_t fullSteps = totalSteps - remainderSteps;

    vector<array<double, 4>> vDiagnostics;
    vDiagnostics.reserve(_params.iGetNumericalSteps() + 1);
    vector<vector<double>> vSpectrumHistory;
    vSpectrumHistory.reserve(_params.iGetNumericalSteps() + 1);

    SolutionData vStateCurrent = vTargetStart; // set phi_{i} = phi(0)
    SolutionData vStateNext(_params, _paths, InitialState); // phi_{i+1}

    SolutionData vNonlinearTermPrevious(_params, _paths, InitialState); // N(phi_{i-1})
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        vNonlinearTermPrevious[i] = complex<double>{0.0, 0.0};

    SolutionData vNonlinearTermCurrent(_params, _paths, InitialState); // N(phi_{i})
    SolutionData vForwardSpatialDerivative1(_params, _paths, InitialState); // (phi_{i})_{x_1}
    SolutionData vForwardSpatialDerivative2(_params, _paths, InitialState); // (phi_{i})_{x_2}

    vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2() });
    vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
    for (size_t i = 1; i < totalSteps + 1; ++i) {
        dTimePoint = i * _params.dTimeStep;
        for (size_t j = 0; j < 4; ++j) {
            
            // phi_{x_1} and phi_{x_2} in Fourier space
            for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                vForwardSpatialDerivative1[k] = _params.vDifferentialOperator1[k] * vStateCurrent[k];
                vForwardSpatialDerivative2[k] = _params.vDifferentialOperator2[k] * vStateCurrent[k];
            }
            
            // phi_{x_1}^2 and phi_{x_2}^2 in physical space
            _fftwPlan.ifft2InPlace(_params, vForwardSpatialDerivative1.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vForwardSpatialDerivative2.getDataPointer());
            checkCFL(vForwardSpatialDerivative1, vForwardSpatialDerivative2);

            // build nonlinear term in physical space
            for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                const double phidx1 = vForwardSpatialDerivative1[k].real();
                const double phidx2 = vForwardSpatialDerivative2[k].real();
                vNonlinearTermCurrent[k] = complex<double>{0.5 * (phidx1 * phidx1 + phidx2 * phidx2), 0.0};
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

        saveForwardState(dTimePoint, i, fullSteps, stepsPerFile, totalSteps, vHistoryIntermediate, vHistoryRemainder, vStateCurrent);
        vDiagnostics.push_back({dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2()});
        vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
    }
    vTargetEnd.moveDataFrom(vStateCurrent);
    setSolutionState(SolveTerminalState, vTargetEnd);
    if (_params.bOptimizeSolution == 0) {
        saveSolutionDiagnostics(vDiagnostics); 
        saveSolutionSpectrum(vSpectrumHistory);
    }
}

void Solver::saveForwardState(double dTimePoint, size_t forwardIndex, size_t fullSteps, size_t stepsPerFile,  size_t totalSteps,
                              SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vStateCurrent) {
    if (forwardIndex <= fullSteps) {
        const size_t localIndex = (forwardIndex - 1) % stepsPerFile;
        vHistoryIntermediate.setData(vStateCurrent, localIndex);
        if (forwardIndex % stepsPerFile == 0)
            vHistoryIntermediate.saveData(IntermediateHistory, dTimePoint);
    }
    else {
        const size_t localIndex = forwardIndex - fullSteps - 1;
        vHistoryRemainder.setData(vStateCurrent, localIndex);
        if (forwardIndex == totalSteps)
            vHistoryRemainder.saveData(RemainderHistory, dTimePoint);
    }
}

void Solver::loadForwardState(double dTimePoint, size_t forwardIndex, size_t fullSteps, size_t stepsPerFile, 
                              SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vForwardStateCurrent) {
    if (forwardIndex > fullSteps) {
        const size_t localIndex = forwardIndex - fullSteps - 1;
        vHistoryRemainder.getData(vForwardStateCurrent, localIndex);
    }
    else {
        const size_t fileIndex  = (forwardIndex - 1) / stepsPerFile;
        const size_t localIndex = (forwardIndex - 1) % stepsPerFile;

        // Reload only when current time point reaches current file range
        if ( forwardIndex == (fileIndex + 1) * stepsPerFile ) 
            vHistoryIntermediate.loadData(IntermediateHistory, dTimePoint);
        vHistoryIntermediate.getData(vForwardStateCurrent, localIndex);
    }
}

void Solver::solveBackwardInTime(SolutionData& vObjectiveGradient, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, const SolutionData& vTargetEnd) {
    double dTimePoint = 0.0;
    const size_t totalSteps = _params.iGetNumericalSteps();
    const size_t stepsPerFile = _params.iGetNumericalStepsPerFile();
    const size_t remainderSteps = totalSteps % stepsPerFile;
    const size_t fullSteps = totalSteps - remainderSteps;

    // set phi^*_{i} = 2 * phi(T)
    SolutionData vStateCurrent = vTargetEnd;
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        vStateCurrent[i] *= 2.0;

    SolutionData vForwardStateCurrent = vTargetEnd; // phi_{i}
    SolutionData vStateNext(_params, _paths, InitialState); // phi^*_{i-1}
    SolutionData vNonlinearTermPrevious(_params, _paths, InitialState); // N(phi^*_{i+1})
    for (size_t i = 0; i < _params.iTotalGridSize; ++i)
        vNonlinearTermPrevious[i] = complex<double>{0.0, 0.0};

    SolutionData vNonlinearTermCurrent(_params, _paths, InitialState); // N(phi^*_{i})
    SolutionData vForwardSpatialDerivative1(_params, _paths, InitialState); // (phi_{i})_{x_1}
    SolutionData vForwardSpatialDerivative2(_params, _paths, InitialState); // (phi_{i})_{x_2}
    SolutionData vForwardLaplacian(_params, _paths, InitialState); // Lap(phi_{i})
    SolutionData vBackwardSpatialDerivative1(_params, _paths, InitialState); // (phi^*_{i})_{x_1}
    SolutionData vBackwardSpatialDerivative2(_params, _paths, InitialState); // (phi^*_{i})_{x_2}

    for (size_t i = totalSteps; i > 0; --i) {
        dTimePoint = i * _params.dTimeStep;
        loadForwardState(dTimePoint, i, fullSteps, stepsPerFile, vHistoryIntermediate, vHistoryRemainder, vForwardStateCurrent);

        // Laplacian(phi), phi_{x_1} and phi_{x_2} in physical space
        for (size_t j = 0; j < _params.iTotalGridSize; ++j) {
            vForwardSpatialDerivative1[j] = _params.vDifferentialOperator1[j] * vForwardStateCurrent[j];
            vForwardSpatialDerivative2[j] = _params.vDifferentialOperator2[j] * vForwardStateCurrent[j];
            vForwardLaplacian[j] = _params.vLaplaceOperator[j] * vForwardStateCurrent[j];
        }
        _fftwPlan.ifft2InPlace(_params,vForwardSpatialDerivative1.getDataPointer());
        _fftwPlan.ifft2InPlace(_params,vForwardSpatialDerivative2.getDataPointer());
        _fftwPlan.ifft2InPlace(_params,vForwardLaplacian.getDataPointer());

        for (size_t j = 0; j < 4; ++j) {
            
            // phi^*_{x_1} and phi^*_{x_2} in Fourier space
            for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                vBackwardSpatialDerivative1[k] = _params.vDifferentialOperator1[k] * vStateCurrent[k];
                vBackwardSpatialDerivative2[k] = _params.vDifferentialOperator2[k] * vStateCurrent[k];
            }
            
            // phi^*, phi^*_{x_1} and phi^*_{x_2} in physical space
            _fftwPlan.ifft2InPlace(_params, vStateCurrent.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vBackwardSpatialDerivative1.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vBackwardSpatialDerivative2.getDataPointer());

            // build nonlinear term in physical space
            for (size_t k = 0; k < _params.iTotalGridSize; ++k) {
                const double adj = vStateCurrent[k].real();
                const double Lap = vForwardLaplacian[k].real();
                const double adjdx1 = vBackwardSpatialDerivative1[k].real();
                const double adjdx2 = vBackwardSpatialDerivative2[k].real();
                const double phidx1 = vForwardSpatialDerivative1[k].real();
                const double phidx2 = vForwardSpatialDerivative2[k].real();
                vNonlinearTermCurrent[k] = complex<double>{-1.0 * (adjdx1 * phidx1 + adjdx2 * phidx2) - (Lap * adj), 0.0};
            }

            // dealias transformed nonlinear term and IMEX RK4 update
            _fftwPlan.fft2InPlace(vStateCurrent.getDataPointer());
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
    }
    vObjectiveGradient.moveDataFrom(vStateCurrent);
    setSolutionState(SolveBackwardInitialState, vObjectiveGradient);
}

void Solver::solveRiemmanianOptimization(double& dObjectiveValue, SolutionData& vObjectiveGradient, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                     SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    
    // Finds the maximum L2 energy via Riemmanian Conjugate Gradient with Polak-Ribiere Momentum
    if (_params.bOptimizeSolution == 1) {
        
        // Clear old data
        filesystem::remove(_paths.fOptimizationDiagnostics);
        filesystem::remove(_paths.fOptimizationLineSearch);
        
        size_t iter = 1;
        size_t iMaxIter = 1000;

        dObjectiveValue = vTargetEnd.getEnergyL2(); 
        double dObjectiveDelta = 0.0; 
        double dStepSize = 0.0; 
        double dManifoldSize = vTargetStart.getEnergyL2(); 
        double dObjectiveGradientSize = vObjectiveGradient.getEnergyL2(); 
        double dMomentumSize = 0.0; 
        array<double, 7> vDiagnostics;
        array<double, 7> vOptimalEnergySolution; 

        vDiagnostics = {dObjectiveValue, dObjectiveDelta, dStepSize, dManifoldSize, _timer.elapsedSeconds(), dObjectiveGradientSize, dMomentumSize};
        saveOptimizationDiagnostics(vDiagnostics);

        cout << left << setw(5)  << "Iter" << setw(18) << "Objective (J)" << setw(10) << "Adjoint" << setw(10) 
             << "Stepsize" << setw(6) << "Brent" << setw(10) << "Forward" << setw(10) << "Objective Delta" << '\n' << string(75, '-') << '\n';

        size_t iMomentumCounter = 0;
        SolutionData RawUpdate(_params, _paths, InitialState); 
        SolutionData DirectionCurr(_params, _paths, InitialState); 
        SolutionData DirectionPrev(_params, _paths, InitialState); 
        SolutionData VectorTransport(_params, _paths, InitialState);
        SolutionData ProjectedObjGradPrev(_params, _paths, InitialState); 
        SolutionData ProjectedObjGradCurr(_params, _paths, InitialState); 
        SolutionData TransportedProjectedObjGradPrev(_params, _paths, InitialState); 
        SolutionData DeltaProjectedObjGrad(_params, _paths, InitialState); 
        double dAngleFwdICAndGradJ = vTargetStart.getInnerProductL2With(vObjectiveGradient);
        double dRawUpdateSize = 0.0;
        double dAngleRawUpdateAndDirectionPrev = 0.0; 
        double dAngleRawUpdateAndProjectedObjGradPrev = 0.0; 
        double dProjectedObjGradPrevSize = 0.0; 
        double dProjectedObjGradCurrSize = 0.0;
        double dAngleProjGradJAndDeltaProjGradJ = 0.0;
        double dDirectionSize = 0.0; 
        double dAscentSize = 0.0;
        double dRetractionSize = 0.0;
        double dObjectiveValuePrev = 0.0;

        dObjectiveDelta = 1.0;
        while ( abs(dObjectiveDelta) > _params.dOptimizationTolerance && iter <= iMaxIter ) {
            cout << left << setw(5)  << iter << setw(18) << setprecision(12) << dObjectiveValue << flush;
            cout << defaultfloat << setprecision(6);
            
            // Solve adjoint equation
            setSolutionInTime(SolveBackwardInTime, vObjectiveGradient, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            _timer.printIterationInterval(); 
            
            // Update initial data via RCG
            dObjectiveGradientSize = vObjectiveGradient.getEnergyL2(); 
            dAngleFwdICAndGradJ = vTargetStart.getInnerProductL2With(vObjectiveGradient);
            ProjectedObjGradCurr = vObjectiveGradient - (dAngleFwdICAndGradJ / _params.dInitialEnergy) * vTargetStart;
            if (iter > 1) {
                dRawUpdateSize = RawUpdate.getEnergyL2();
                dAngleRawUpdateAndDirectionPrev = RawUpdate.getInnerProductL2With(DirectionPrev);                       
                dAngleRawUpdateAndProjectedObjGradPrev = RawUpdate.getInnerProductL2With(ProjectedObjGradPrev);          
                VectorTransport = (DirectionPrev - (dAngleRawUpdateAndDirectionPrev / dRawUpdateSize) * (RawUpdate)) 
                                  * sqrt(_params.dInitialEnergy / dRawUpdateSize);          
                TransportedProjectedObjGradPrev = (ProjectedObjGradPrev - (dAngleRawUpdateAndProjectedObjGradPrev / dRawUpdateSize) * (RawUpdate))
                                  * sqrt(dManifoldSize / dRawUpdateSize); 
                DeltaProjectedObjGrad = ProjectedObjGradCurr - TransportedProjectedObjGradPrev;                       
                dAngleProjGradJAndDeltaProjGradJ = ProjectedObjGradCurr.getInnerProductL2With(DeltaProjectedObjGrad);
                dProjectedObjGradPrevSize = ProjectedObjGradPrev.getEnergyL2();     
                dMomentumSize = dAngleProjGradJAndDeltaProjGradJ / dProjectedObjGradPrevSize;                    
                iMomentumCounter = iMomentumCounter + 1;
                if (iter % 20 == 0)
                    dMomentumSize = 0.0;
                DirectionCurr = ProjectedObjGradCurr + (dMomentumSize * VectorTransport); 
            }
            else {
                dStepSize = 1e5;
                DirectionCurr = ProjectedObjGradCurr;
            }
                         
            dDirectionSize = DirectionCurr.getNormL2();                       
            dProjectedObjGradCurrSize = ProjectedObjGradCurr.getNormL2();
            dAscentSize = DirectionCurr.getInnerProductL2With(ProjectedObjGradCurr) / ( dDirectionSize * dProjectedObjGradCurrSize);
            if (dAscentSize < 0)
                DirectionCurr = ProjectedObjGradCurr;             

            solveLineSearchOptimization(dStepSize, DirectionCurr, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            if (dStepSize == 0)
                dStepSize = dAngleFwdICAndGradJ / dObjectiveGradientSize; 
                                     
            RawUpdate = vTargetStart + ( dStepSize * DirectionCurr );                                   
            dRetractionSize =  sqrt(_params.dInitialEnergy) / RawUpdate.getNormL2();                                                
            vTargetStart = dRetractionSize * RawUpdate;                                                            
            dManifoldSize = vTargetStart.getEnergyL2();                                   

            
            // Solve updated forward problem
            setSolutionInTime(SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            vTargetStart.saveData(InitialState);
            _timer.printIterationInterval(); 

            dObjectiveValuePrev = dObjectiveValue;                                                                                
            dObjectiveValue = vTargetEnd.getEnergyL2();                                                  
            iter = iter + 1;                                                                               
            ProjectedObjGradPrev = ProjectedObjGradCurr;                                                                  
            DirectionPrev = DirectionCurr;                                                                            
                                              
            dObjectiveDelta = (dObjectiveValue - dObjectiveValuePrev)/(dObjectiveValuePrev);                                                                                                                                                     

            cout << setw(10) << dObjectiveDelta << flush << '\n';
            vDiagnostics = {dObjectiveValue, dObjectiveDelta, dStepSize, dManifoldSize, _timer.elapsedSeconds(), dObjectiveGradientSize, dMomentumSize};
            saveOptimizationDiagnostics(vDiagnostics);
        }

        _params.bOptimizeSolution = 0;
        setSolutionInTime(SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);

        EnergyData maxEnergyResult = getMaxEnergyL2InTimeWindow();
        vOptimalEnergySolution = { _params.dInitialEnergy, _params.dDomainFactor1, _params.dDomainFactor2, _params.dTimeWindow, dObjectiveValue, 
                                   maxEnergyResult.dTimepoint, maxEnergyResult.dEnergy };
        saveSolutionBranchAndPowerLaws(vOptimalEnergySolution);
        cout << string(75, '-') << '\n';
        _timer.printInterval("Energy maximization problem solved at ");
    }
}

void Solver::solveLineSearchOptimization(double& dStepSize, SolutionData& DirectionCurr, SolutionData& vTargetStart, 
                                         SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    
    // Finds the optimal step size using bracketing followed by Brent's method.
    constexpr double relativeTolerance    = 1e-5;
    constexpr double goldenExpansionFactor = 1.618034;
    constexpr double maximumExpansionFactor = 100.0;
    constexpr double goldenSectionFactor  = 0.381966;
    constexpr double numericalEpsilon   = 1e-10;

    size_t iteration = 0;
    constexpr size_t maximumIterations = 1000;

    vector<double> vLineSearchHistory;
    vLineSearchHistory.reserve(maximumIterations);

    SolutionData RawUpdate(_params, _paths, InitialState);
    SolutionData RetractedState(_params, _paths, InitialState);
    
    // Evaluate f(step) = -J(step), where J is the terminal L2 energy. Therefore, minimizing f maximizes J.
    auto evaluateStep = [&](double step) {
        RawUpdate = vTargetStart + step * DirectionCurr;
        const double rawNorm = RawUpdate.getNormL2();
        if (rawNorm == 0.0 || !isfinite(rawNorm))
            return numeric_limits<double>::infinity();
        const double retraction = sqrt(_params.dInitialEnergy) / rawNorm;
        RetractedState = retraction * RawUpdate;
        setSolutionInTime( SolveForwardInTime, RetractedState, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
        const double objectiveValue = vTargetEnd.getEnergyL2();
        if (iteration < maximumIterations) {
            vLineSearchHistory.push_back(objectiveValue);
            ++iteration;
        }
        return -objectiveValue;
    };

    // Bracketing stage
    double previousStepSize = 0.0;
    double middleStepSize = dStepSize;
    double previousObjective = evaluateStep(previousStepSize);
    double middleObjective = evaluateStep(middleStepSize);

    if (middleObjective > previousObjective) {
        swap(previousStepSize, middleStepSize);
        swap(previousObjective, middleObjective);
    }

    double nextStepSize = middleStepSize + goldenExpansionFactor * (middleStepSize - previousStepSize);
    double nextObjective = evaluateStep(nextStepSize);

    while (middleObjective >= nextObjective && iteration < maximumIterations) {
        const double interpolationTermPrevious = (middleStepSize - previousStepSize) * (middleObjective - nextObjective);
        const double interpolationTermNext = (middleStepSize - nextStepSize) * (middleObjective - previousObjective);
        const double interpolationDifference = interpolationTermNext - interpolationTermPrevious;
        const double interpolationDenominator = 2.0 * copysign(max(abs(interpolationDifference), numericalEpsilon), interpolationDifference);
        double trialStepSize = middleStepSize - ((middleStepSize - nextStepSize) * interpolationTermNext - (middleStepSize - previousStepSize) * interpolationTermPrevious) / interpolationDenominator;
        const double maximumTrialStepSize = middleStepSize + maximumExpansionFactor * (nextStepSize - middleStepSize);
        double trialObjective = 0.0;

        if ((middleStepSize - trialStepSize) * (trialStepSize - nextStepSize) > 0.0) {
            trialObjective = evaluateStep(trialStepSize);

            if (trialObjective < nextObjective) {
                previousStepSize = middleStepSize;
                previousObjective = middleObjective;
                middleStepSize = trialStepSize;
                middleObjective = trialObjective;
                break;
            }

            if (trialObjective > middleObjective) {
                nextStepSize = trialStepSize;
                nextObjective = trialObjective;
                break;
            }

            trialStepSize = nextStepSize + goldenExpansionFactor * (nextStepSize - middleStepSize);
            trialObjective = evaluateStep(trialStepSize);
        }
        else if ((nextStepSize - trialStepSize) * (trialStepSize - maximumTrialStepSize) > 0.0) {
            trialObjective = evaluateStep(trialStepSize);

            if (trialObjective < nextObjective) {
                middleStepSize = nextStepSize;
                middleObjective = nextObjective;
                nextStepSize = trialStepSize;
                nextObjective = trialObjective;
                trialStepSize = nextStepSize + goldenExpansionFactor * (nextStepSize - middleStepSize);
                trialObjective = evaluateStep(trialStepSize);
            }
        }
        else if ((trialStepSize - maximumTrialStepSize) * (maximumTrialStepSize - nextStepSize) >= 0.0) {
            trialStepSize = maximumTrialStepSize;
            trialObjective = evaluateStep(trialStepSize);
        }
        else {
            trialStepSize = nextStepSize + goldenExpansionFactor * (nextStepSize - middleStepSize);
            trialObjective = evaluateStep(trialStepSize);
        }

        previousStepSize = middleStepSize;
        middleStepSize = nextStepSize;
        nextStepSize = trialStepSize;
        previousObjective = middleObjective;
        middleObjective = nextObjective;
        nextObjective = trialObjective;
    }

    if (iteration >= maximumIterations) {
        saveLineSearch(vLineSearchHistory);
        _timer.printIterationInterval();
        cout << setw(6) << iteration << flush;
        return;
    }

    // Brent minimization stage
    double bracketLowerBound = min(previousStepSize, nextStepSize);
    double bracketUpperBound = max(previousStepSize, nextStepSize);
    double secondPreviousBestStep = middleStepSize;
    double previousBestStep = middleStepSize;
    double bestStep = middleStepSize;
    double bestObjective = middleObjective;
    double previousBestObjective = middleObjective;
    double secondPreviousBestObjective = middleObjective;
    double proposedStepOffset = 0.0;
    double previousStepOffset = 0.0;

    while (iteration < maximumIterations) {
        const double bracketMidpoint = 0.5 * (bracketLowerBound + bracketUpperBound);
        const double absoluteTolerance = relativeTolerance * abs(bestStep) + numericalEpsilon;
        const double doubledTolerance = 2.0 * absoluteTolerance;

        if (abs(bestStep - bracketMidpoint) <= doubledTolerance - 0.5 * (bracketUpperBound - bracketLowerBound))
            break;

        bool useParabolicStep = false;
        if (abs(previousStepOffset) > absoluteTolerance) {
            const double interpolationTermPrevious = (bestStep - previousBestStep) * (bestObjective - secondPreviousBestObjective);
            double interpolationTermSecondPrevious = (bestStep - secondPreviousBestStep) * (bestObjective - previousBestObjective);
            double parabolicNumerator = (bestStep - secondPreviousBestStep) * interpolationTermSecondPrevious - (bestStep - previousBestStep) * interpolationTermPrevious;
            double parabolicDenominator = 2.0 * (interpolationTermSecondPrevious - interpolationTermPrevious);

            if (parabolicDenominator > 0.0)
                parabolicNumerator = -parabolicNumerator;

            parabolicDenominator = abs(parabolicDenominator);
            const double savedPreviousStepOffset = previousStepOffset;
            previousStepOffset = proposedStepOffset;

            if (parabolicDenominator > 0.0 && abs(parabolicNumerator) < abs(0.5 * parabolicDenominator * savedPreviousStepOffset) && parabolicNumerator > parabolicDenominator * (bracketLowerBound - bestStep) && parabolicNumerator < parabolicDenominator * (bracketUpperBound - bestStep)) {
                proposedStepOffset = parabolicNumerator / parabolicDenominator;
                useParabolicStep = true;

                const double proposedTrialStep = bestStep + proposedStepOffset;

                if (proposedTrialStep - bracketLowerBound < doubledTolerance || bracketUpperBound - proposedTrialStep < doubledTolerance)
                    proposedStepOffset = copysign(absoluteTolerance, bracketMidpoint - bestStep);
            }
        }

        if (!useParabolicStep) {
            previousStepOffset = (bestStep >= bracketMidpoint) ? bracketLowerBound - bestStep : bracketUpperBound - bestStep;
            proposedStepOffset = goldenSectionFactor * previousStepOffset;
        }

        const double trialStep = abs(proposedStepOffset) >= absoluteTolerance ? bestStep + proposedStepOffset : bestStep + copysign(absoluteTolerance, proposedStepOffset);
        const double trialObjective = evaluateStep(trialStep);

        if (trialObjective <= bestObjective) {
            if (trialStep >= bestStep)
                bracketLowerBound = bestStep;
            else
                bracketUpperBound = bestStep;

            secondPreviousBestStep = previousBestStep;
            secondPreviousBestObjective = previousBestObjective;
            previousBestStep = bestStep;
            previousBestObjective = bestObjective;
            bestStep = trialStep;
            bestObjective = trialObjective;
        }
        else {
            if (trialStep < bestStep)
                bracketLowerBound = trialStep;
            else
                bracketUpperBound = trialStep;

            if (trialObjective <= previousBestObjective || previousBestStep == bestStep) {
                secondPreviousBestStep = previousBestStep;
                secondPreviousBestObjective = previousBestObjective;
                previousBestStep = trialStep;
                previousBestObjective = trialObjective;
            }
            else if (trialObjective <= secondPreviousBestObjective || secondPreviousBestStep == bestStep || secondPreviousBestStep == previousBestStep) {
                secondPreviousBestStep = trialStep;
                secondPreviousBestObjective = trialObjective;
            }
        }
    }

    dStepSize = isfinite(bestStep) ? bestStep : 0.0;
    saveLineSearch(vLineSearchHistory);
    _timer.printIterationInterval();
    cout << setw(6) << iteration << flush;
}

// Public functions
Solver::Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer) 
                : _params(params), _paths(paths), _fftwPlan(fftwPlan), _timer(timer) {};

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
                findContinuationForInitialData(vTargetState);
            }
            vTargetState.saveData(InitialState);
            break;
        }

        case SolveTerminalState:
            vTargetState.saveData(TerminalState);
            break;
        
        case SolveBackwardInitialState: 
            vTargetState.saveData(BackwardInitialState);
            break;
    }  
}

void Solver::setSolutionInTime(TimeSolutionType targetType, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        
        case SolveForwardInTime: 
            solveForwardInTime(vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            break;

        case SolveBackwardInTime:
            solveBackwardInTime(vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            break;
    }  
}

double Solver::getOptimalSolution(OptimizeSolutionType targetType, SolutionData& vObjectiveGradient, SolutionData& vTargetStart, 
                SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {

    // return energy or step-size for respective cases
    double dTargetValue = 0.0; 
    switch (targetType) {
        
        case OptimizeEnergyAmplification:  
            cout << "\n";    
            solveRiemmanianOptimization(dTargetValue, vObjectiveGradient, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            cout << "\n";
            break;

        case OptimizeLineSearchStepSize:
            solveLineSearchOptimization(dTargetValue, vObjectiveGradient, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            break;

    } 
    return dTargetValue;
}
