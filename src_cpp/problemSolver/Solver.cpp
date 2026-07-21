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

// Private functions
void Solver::saveSolutionDiagnostics(const vector<array<double, 4>>& vDiagnostics) {
    if (!_mpi.isRoot()) {
        return;
    }
    
    ofstream file(_paths.getEnergyEvolutionFile());
    file << setprecision(16) << scientific;

    for (const auto& row : vDiagnostics) {
        file << row[0] << ' ' << row[1] << ' ' << row[2] << ' ' << row[3] << '\n';
    }
}

void Solver::saveSolutionSpectrum(const vector<vector<double>>& vSpectrumHistory) {
    if (!_mpi.isRoot()) {
        return;
    }
    
    ofstream file(_paths.getFourierSpectrumEvolutionFile());
    file << setprecision(16) << scientific;
    const size_t maxRadialBin = _params.getMaxRadialBin();
    const double radialBinWidth = _params.getRadialBinWidth();

    for (size_t i = 0; i <= maxRadialBin; ++i) {
        const double radialBinRadius = static_cast<double>(i) * radialBinWidth;
        file << radialBinRadius;
        
        for (const auto& spectrum : vSpectrumHistory) {
            file << ' ' << spectrum[i];
        }
        file << '\n';
    }
}

void Solver::saveOptimizationDiagnostics(const array<double, 7>& vDiagnostics) {
    if (!_mpi.isRoot()) {
        return;
    }
    
    ofstream file(_paths.getOptimizationDiagnosticsFile(), ios::app);
    file << setprecision(16) << scientific;
    file << vDiagnostics[0] << ' ' << vDiagnostics[1] << ' ' << vDiagnostics[2] << ' ' << vDiagnostics[3] 
         << ' ' << vDiagnostics[4] << ' ' << vDiagnostics[5] << ' ' << vDiagnostics[6] << '\n';
}

void Solver::saveSolutionBranchAndPowerLaws(const array<double, 7>& vOptimalEnergySolution) {
    if (!_mpi.isRoot()) {
        return;
    }
    
    vector<array<double, 7>> vBranch;
    vector<array<double, 7>> vK;
    vector<array<double, 7>> vL;
    vector<array<double, 7>> vTK;
    vector<array<double, 7>> vTL;

    {
        ifstream inputBranch(_paths.getSolutionBranchesFile());
        ifstream inputK(_paths.getInitialEnergyPowerLawFile());
        ifstream inputL(_paths.getDomainSizePowerLawFile());
        ifstream inputTK(_paths.getEnergyTimeWindowPowerLawFile());
        ifstream inputTL(_paths.getDomainTimeWindowPowerLawFile());
        array<double, 7> solutionBranch;
        array<double, 7> solutionK;
        array<double, 7> solutionL;
        array<double, 7> solutionTK;
        array<double, 7> solutionTL;

        while (inputBranch >> solutionBranch[0] >> solutionBranch[1] >> solutionBranch[2] >> solutionBranch[3] >> solutionBranch[4] >> solutionBranch[5] >> solutionBranch[6]) {
            vBranch.push_back(solutionBranch);
        }
        while (inputK >> solutionK[0] >> solutionK[1] >> solutionK[2] >> solutionK[3] >> solutionK[4] >> solutionK[5] >> solutionK[6]) {
            vK.push_back(solutionK);
        }
        while (inputL >> solutionL[0] >> solutionL[1] >> solutionL[2] >> solutionL[3] >> solutionL[4] >> solutionL[5] >> solutionL[6]) {
            vL.push_back(solutionL);
        }
        while (inputTK >> solutionTK[0] >> solutionTK[1] >> solutionTK[2] >> solutionTK[3] >> solutionTK[4] >> solutionTK[5] >> solutionTK[6]) {
            vTK.push_back(solutionTK);
        }
        while (inputTL >> solutionTL[0] >> solutionTL[1] >> solutionTL[2] >> solutionTL[3] >> solutionTL[4] >> solutionTL[5] >> solutionTL[6]) {
            vTL.push_back(solutionTL);
        }
    }

    auto updateAndSave = [&](vector<array<double, 7>>& solutions, const string& filename, auto matches, auto sortingRule) {
        auto existing = find_if(solutions.begin(), solutions.end(), [&](const auto& solution) { return matches(solution, vOptimalEnergySolution);});

        if (existing == solutions.end()) {
            solutions.push_back(vOptimalEnergySolution);
        }
        else if (vOptimalEnergySolution[6] > (*existing)[6]) {
            *existing = vOptimalEnergySolution;
        }

        sort(solutions.begin(), solutions.end(), sortingRule);
        ofstream output(filename);
        output << scientific << setprecision(16);

        for (const auto& solution : solutions) {
            size_t solutionSize = solution.size();
            for (size_t column = 0; column < solutionSize; ++column) {
                if (column > 0) {
                    output << ' ';
                }
                output << solution[column];
            }
            output << '\n';
        }
    };

    updateAndSave(vBranch, _paths.getSolutionBranchesFile(),
        [](const auto& a, const auto& b) {
            return a[3] == b[3];
        },
        [](const auto& a, const auto& b) {
            return a[3] < b[3];
        });

    updateAndSave(vK, _paths.getInitialEnergyPowerLawFile(),
        [](const auto& a, const auto& b) {
            return a[0] == b[0];
        },
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    updateAndSave(vL, _paths.getDomainSizePowerLawFile(),
        [](const auto& a, const auto& b) {
            return a[1] == b[1] && a[2] == b[2];
        },
        [](const auto& a, const auto& b) {
            if (a[1] != b[1])
                return a[1] < b[1];
            return a[2] < b[2];
        });

    updateAndSave(vTK, _paths.getEnergyTimeWindowPowerLawFile(),
        [](const auto& a, const auto& b) {
            return a[0] == b[0];
        },
        [](const auto& a, const auto& b) {
            return a[0] < b[0];
        });

    updateAndSave(vTL, _paths.getDomainTimeWindowPowerLawFile(),
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
    if (!_mpi.isRoot()) {
        return;
    }
    
    ofstream file(_paths.getOptimizationLineSearchFile(), ios::app);
    const double dNaN = numeric_limits<double>::quiet_NaN();
    file << setprecision(16) << scientific;

    for (size_t i = 0; i < 1000; ++i) {
        if (i < vLineSearchHistory.size()) {
            file << vLineSearchHistory[i];
        }
        else {
            file << dNaN;
        }

        if (i + 1 < 1000) {
            file << ' ';
        }
    }
    file << '\n';
}

EnergyData Solver::getMaxEnergyL2InTimeWindow() {
    EnergyData result{-numeric_limits<double>::infinity(), 0.0};

    // Ensure rank 0 has finished writing and closing the file.
    MPI_Barrier(MPI_COMM_WORLD);

    if (_mpi.isRoot()) {
        ifstream inputFile(_paths.getEnergyEvolutionFile());
        double timepoint;
        double energyL2;
        double energyH1;
        double energyH2;

        while (inputFile >> timepoint >> energyL2 >> energyH1 >> energyH2) {
            if (energyL2 > result.dEnergy) {
                result.dEnergy = energyL2;
                result.dTimepoint = timepoint;
            }
        }
    }

    double resultValues[2] = {result.dEnergy, result.dTimepoint};
    MPI_Bcast(resultValues, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    result.dEnergy = resultValues[0];
    result.dTimepoint = resultValues[1];

    return result;
}

void Solver::checkCFL(SolutionData& vData1, SolutionData& vData2) {
    constexpr double CFL_LIMIT = 0.9;
    double localMaxGradientMagnitude = 0.0;
    const size_t localGridSize = _mpi.getLocalGridSize();

    for (size_t i = 0; i < localGridSize; ++i) {
        const double velocity = hypot(vData1[i].real(), vData2[i].real()); // sqrt(ux * ux + uy * uy)
        localMaxGradientMagnitude = max(localMaxGradientMagnitude, velocity);
    }

    double globalMaxGradientMagnitude = 0.0;
    MPI_Allreduce(&localMaxGradientMagnitude, &globalMaxGradientMagnitude, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    const double cflValue = _params.getTimeStep() * globalMaxGradientMagnitude / _params.getSpaceStep();

    if (cflValue > CFL_LIMIT) {
        if (_mpi.isRoot()) {
            cerr << "ERROR: CFL condition violated." << " Maximum gradient magnitude = " << globalMaxGradientMagnitude 
                << ", CFL value = " << cflValue << ", CFL Limit = " << CFL_LIMIT << '\n';
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
}

void Solver::setInitialCondition(SolutionData& vTargetState) {  
    const size_t gridSize1 = _params.getGridSize1();
    const size_t localGridSize2 = _mpi.getLocalGridSize2();
    const size_t localGridStart2 = _mpi.getLocalGridStart2();

    const vector<double>& grid1 = _params.getGrid1();
    const vector<double>& grid2 = _params.getGrid2();
    const double domainFactor1 = _params.getDomainFactor1();
    const double domainFactor2 = _params.getDomainFactor2();

    if (_params.getInitialGuessName() == "s1") {
        size_t stateNumber = 0;

        for (size_t i = 0; i < localGridSize2; ++i) {
            for (size_t j = 0; j < gridSize1; ++j) {
                const size_t k = _params.getIndex(i, j, gridSize1);
                const double x = grid1[k];
                const double y = grid2[k];
                vTargetState(j, i, stateNumber) = complex<double>( sin(x / domainFactor1 + y / domainFactor2)
                                                                 + sin(x / domainFactor1) + sin(y / domainFactor2), 0.0);
            }
        }
        _fftwPlan.fft2InPlace(vTargetState.getDataPointer());
    }
    else if (_params.getInitialGuessName() == "randfour") {
        const size_t stateNumber = 0;
        const double kmax  = max(_params.getGridSize1(), _params.getGridSize2()) / 2.0;
        const double delta = 6.0 / kmax;

        mt19937 rng(random_device{}());
        uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);

        complex<double>* vTempState = vTargetState.getStatePointer(stateNumber);
        const vector<double>& modes1 = _params.getWavenumbersNonlinear1();
        const vector<double>& modes2 = _params.getWavenumbersNonlinear2();

        // Build randomized Fourier coefficients
        for (size_t iLocal = 0; iLocal < localGridSize2; ++iLocal) {
            const size_t i = localGridStart2 + iLocal;

            for (size_t j = 0; j < gridSize1; ++j) {
                const size_t k = _params.getIndex(iLocal, j, gridSize1);
                const double k1 = modes1[j];
                const double k2 = modes2[i];
                const double kabs = sqrt(k1 * k1 + k2 * k2);
                const double A = exp(-delta * kabs);
                const double phi = dist(rng);
                vTempState[k] = A * complex<double>(cos(phi), sin(phi));
            }
        }

        // zero mean mode
        if (localGridStart2 == 0 && localGridSize2 > 0) {
            vTempState[0] = complex<double>(0.0, 0.0);
        }
        
        _fftwPlan.ifft2InPlace(_params, vTempState);  // real(ifft2(uhat))

        const size_t localGridSize = _mpi.getLocalGridSize();
        for (size_t i = 0; i < localGridSize; ++i) {
            vTempState[i] = complex<double>(vTempState[i].real(), 0.0);
        }
        _fftwPlan.fft2InPlace(vTempState);
    }
    else {
        if (_mpi.isRoot()) {
            cerr << "ERROR: Unknown initial condition: " << _params.getInitialGuessName() << endl;
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    vTargetState.setInitialEnergyL2();
}

void Solver::findContinuationForInitialData(SolutionData& vTargetState) {
    filesystem::path fContinuedFile;
    filesystem::path fTempFile;
    const filesystem::path dirInitialData = _paths.getDirInitialData();
    const string testcaseGenericTime = _paths.getTestcaseGenericTime();
    const double optimalTimeWindow = _params.getOptimalTimeWindow();
    double dBestT = _params.getOptimalTimeWindow() + 1.0;            

    // Search for nearest previous time window in directory
    for (const auto& entry : filesystem::directory_iterator(dirInitialData)) {
        
        if (!entry.is_regular_file()) {
            continue;
        }
        
        string filename = entry.path().filename().string();
        if (filename.find(testcaseGenericTime) == string::npos) {
            continue;
        }
        
        size_t pos1 = filename.find("_T_"); 
        if (pos1 == string::npos) {
            continue;
        }
        
        pos1 += 3; // skip "_T_"
        size_t pos2 = filename.find("_opt_", pos1); 
        if (pos2 == string::npos) {
            continue;
        }
        
        double T = stod(filename.substr(pos1, pos2 - pos1)); // Extract T
        if (T <= optimalTimeWindow && ( T > dBestT || dBestT > optimalTimeWindow)) {
            dBestT = T;
            fContinuedFile = entry.path();
        }
    }

    if (dBestT <= optimalTimeWindow) {
        fTempFile = _paths.getInitialDataFile();
        _paths.setInitialDataFile(fContinuedFile);
        vTargetState.loadData(InitialState);
        if (_mpi.isRoot()) {
            cout << "Loaded continued initial data from T = " << dBestT << "\n" << flush;
        }
        _paths.setInitialDataFile(fTempFile); 
    }
    else {
        setInitialCondition(vTargetState);
        if (_mpi.isRoot()) {
            cout << "No saved data, generated initial guess.\n" << flush;
        }
    }
}

void Solver::solveForwardInTime(SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    double dTimePoint = 0.0;
    const size_t totalSteps = _params.getNumericalSteps();
    const size_t stepsPerFile = _params.getNumericalStepsPerFile();
    const size_t remainderSteps = totalSteps % stepsPerFile;
    const size_t fullSteps = totalSteps - remainderSteps;
    const size_t savedStateCount = min(_params.getSavedStates(), totalSteps + 1);
    const size_t savedStepsPerFile = min(stepsPerFile, savedStateCount);
    const size_t savedRemainderSteps = savedStateCount % savedStepsPerFile;
    const size_t savedFullSteps = savedStateCount - savedRemainderSteps;
    size_t savedStateIndex = 0;
    size_t nextSavedStep = 0;    

    const size_t localGridSize = _mpi.getLocalGridSize();
    const size_t localGridSize2 = _mpi.getLocalGridSize2();
    const size_t localGridStart2 = _mpi.getLocalGridStart2();
    const double dt = _params.getTimeStep();
    const vector<complex<double>>& differentialOperator1 = _params.getDifferentialOperator1();
    const vector<complex<double>>& differentialOperator2 = _params.getDifferentialOperator2();
    const vector<complex<double>>& linearOperator = _params.getLinearOperator();
    const vector<complex<double>>& dealiasingOperator = _params.getDealiasingOperator();
    const array<double, 4>& coeffAlphaI = _params.getCoeffAlphaI();
    const array<double, 4>& coeffBetaI = _params.getCoeffBetaI();
    const array<double, 4>& coeffAlphaE = _params.getCoeffAlphaE();
    const array<double, 4>& coeffBetaE = _params.getCoeffBetaE();

    const bool optimizeSolution = _params.getOptimizeSolution();
    const bool activeLineSearch = _params.getActiveLineSearch();

    vector<array<double, 4>> vDiagnostics;
    vDiagnostics.reserve(_params.getNumericalSteps() + 1);
    vector<vector<double>> vSpectrumHistory;
    vSpectrumHistory.reserve(savedStateCount + 1);

    SolutionData vStateCurrent = vTargetStart; // set phi_{i} = phi(0)
    SolutionData vStateNext(_params, _paths, _mpi, InitialState); // phi_{i+1}
    SolutionData vNonlinearTermCurrent(_params, _paths, _mpi, InitialState); // N(phi_{i})
    SolutionData vForwardSpatialDerivative1(_params, _paths, _mpi, InitialState); // (phi_{i})_{x_1}
    SolutionData vForwardSpatialDerivative2(_params, _paths, _mpi, InitialState); // (phi_{i})_{x_2}
    SolutionData vNonlinearTermPrevious(_params, _paths, _mpi, InitialState); // N(phi_{i-1})
    for (size_t i = 0; i < localGridSize; ++i) {
        vNonlinearTermPrevious[i] = complex<double>{0.0, 0.0};
    }

    if (!optimizeSolution && savedStateCount > 0) {
        ++savedStateIndex;
        saveForwardState( dTimePoint, savedStateIndex, savedFullSteps, savedStepsPerFile, savedStateCount, vHistoryIntermediate, vHistoryRemainder, vStateCurrent);
        if (savedStateIndex < savedStateCount) {
            nextSavedStep = static_cast<size_t>(std::llround(static_cast<double>(savedStateIndex) * totalSteps / (savedStateCount - 1)));
        }

        vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2() });
        vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
    }
    else if (!optimizeSolution) {
        vDiagnostics.push_back({ dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2() });
    }
    
    for (size_t i = 1; i < totalSteps + 1; ++i) {
        dTimePoint = i * dt;
        for (size_t j = 0; j < 4; ++j) {
            
            // phi_{x_1} and phi_{x_2} in Fourier space
            for (size_t k = 0; k < localGridSize; ++k) {
                vForwardSpatialDerivative1[k] = differentialOperator1[k] * vStateCurrent[k];
                vForwardSpatialDerivative2[k] = differentialOperator2[k] * vStateCurrent[k];
            }
            
            // phi_{x_1}^2 and phi_{x_2}^2 in physical space
            _fftwPlan.ifft2InPlace(_params, vForwardSpatialDerivative1.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vForwardSpatialDerivative2.getDataPointer());
            checkCFL(vForwardSpatialDerivative1, vForwardSpatialDerivative2);

            // build nonlinear term in physical space
            for (size_t k = 0; k < localGridSize; ++k) {
                const double phidx1 = vForwardSpatialDerivative1[k].real();
                const double phidx2 = vForwardSpatialDerivative2[k].real();
                vNonlinearTermCurrent[k] = complex<double>{0.5 * (phidx1 * phidx1 + phidx2 * phidx2), 0.0};
            }

            // dealias transformed nonlinear term and IMEX RK4 update
            _fftwPlan.fft2InPlace(vNonlinearTermCurrent.getDataPointer());
            for (size_t k = 0; k < localGridSize; ++k) {
                const complex<double> Lin = linearOperator[k];
                vNonlinearTermCurrent[k] *= dealiasingOperator[k];
                vStateNext[k] = ( (1.0 - dt * coeffBetaI[j] * Lin) * vStateCurrent[k]
                                 - dt * coeffAlphaE[j] * vNonlinearTermCurrent[k]
                                 - dt * coeffBetaE[j]  * vNonlinearTermPrevious[k] ) 
                                 / (1.0 + dt * coeffAlphaI[j] * Lin);
            }

            vStateCurrent.swapDataFrom(vStateNext);
            vNonlinearTermPrevious.swapDataFrom(vNonlinearTermCurrent);
        }
        
        // mean-zero correction
        if (localGridStart2 == 0 && localGridSize2 > 0) {
            vStateCurrent[0] = complex<double>{0.0, 0.0};
        }

        if (optimizeSolution && !activeLineSearch) {
            saveForwardState(dTimePoint, i, fullSteps, stepsPerFile, totalSteps, vHistoryIntermediate, vHistoryRemainder, vStateCurrent);
        }
        else if (!optimizeSolution && savedStateIndex < savedStateCount && i == nextSavedStep) {
            ++savedStateIndex;
            saveForwardState(dTimePoint, savedStateIndex, savedFullSteps, savedStepsPerFile, savedStateCount, vHistoryIntermediate, vHistoryRemainder, vStateCurrent);
            
            if (savedStateIndex < savedStateCount) {
                nextSavedStep = static_cast<size_t>(llround(static_cast<double>(savedStateIndex) * totalSteps / (savedStateCount - 1)));
            }

            vDiagnostics.push_back({dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2()});
            vSpectrumHistory.push_back(vStateCurrent.getRadialSpectrum());
        }
        else if (!optimizeSolution) {
            vDiagnostics.push_back({dTimePoint, vStateCurrent.getEnergyL2(), vStateCurrent.getEnergyH1(), vStateCurrent.getEnergyH2()});
        }
    }

    vTargetEnd.moveDataFrom(vStateCurrent);
    setSolutionState(SolveTerminalState, vTargetEnd);
    
    if (!optimizeSolution) {
        saveSolutionDiagnostics(vDiagnostics); 
        saveSolutionSpectrum(vSpectrumHistory);
    }
}

void Solver::saveForwardState(double dTimePoint, size_t forwardIndex, size_t fullSteps, size_t stepsPerFile,  size_t totalSteps,
                              SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vStateCurrent) {
    if (forwardIndex <= fullSteps) {
        const size_t localIndex = (forwardIndex - 1) % stepsPerFile;
        vHistoryIntermediate.setData(vStateCurrent, localIndex);

        if (forwardIndex % stepsPerFile == 0) {
            vHistoryIntermediate.saveData(IntermediateHistory, dTimePoint);
        }
    }
    else {
        const size_t localIndex = forwardIndex - fullSteps - 1;
        vHistoryRemainder.setData(vStateCurrent, localIndex);

        if (forwardIndex == totalSteps) {
            vHistoryRemainder.saveData(RemainderHistory, dTimePoint);
        }
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
        if ( forwardIndex == (fileIndex + 1) * stepsPerFile ) {
            vHistoryIntermediate.loadData(IntermediateHistory, dTimePoint);
        }
        vHistoryIntermediate.getData(vForwardStateCurrent, localIndex);
    }
}

void Solver::solveBackwardInTime(SolutionData& vObjectiveGradient, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, const SolutionData& vTargetEnd) {
    double dTimePoint = 0.0;
    const size_t totalSteps = _params.getNumericalSteps();
    const size_t stepsPerFile = _params.getNumericalStepsPerFile();
    const size_t remainderSteps = totalSteps % stepsPerFile;
    const size_t fullSteps = totalSteps - remainderSteps;

    const size_t localGridSize = _mpi.getLocalGridSize();
    const size_t localGridSize2 = _mpi.getLocalGridSize2();
    const size_t localGridStart2 = _mpi.getLocalGridStart2();
    const double dt = _params.getTimeStep();
    const vector<complex<double>>& differentialOperator1 = _params.getDifferentialOperator1();
    const vector<complex<double>>& differentialOperator2 = _params.getDifferentialOperator2();
    const vector<complex<double>>& linearOperator = _params.getLinearOperator();
    const vector<complex<double>>& laplaceOperator = _params.getLaplaceOperator();
    const vector<complex<double>>& dealiasingOperator = _params.getDealiasingOperator();
    const array<double, 4>& coeffAlphaI = _params.getCoeffAlphaI();
    const array<double, 4>& coeffBetaI = _params.getCoeffBetaI();
    const array<double, 4>& coeffAlphaE = _params.getCoeffAlphaE();
    const array<double, 4>& coeffBetaE = _params.getCoeffBetaE();

    // set phi^*_{i} = 2 * phi(T)
    SolutionData vStateCurrent = vTargetEnd;
    for (size_t i = 0; i < localGridSize; ++i) {
        vStateCurrent[i] *= 2.0;
    }

    SolutionData vForwardStateCurrent = vTargetEnd; // phi_{i}
    SolutionData vStateNext(_params, _paths, _mpi, InitialState); // phi^*_{i-1}
    SolutionData vNonlinearTermCurrent(_params, _paths, _mpi, InitialState); // N(phi^*_{i})
    SolutionData vForwardSpatialDerivative1(_params, _paths, _mpi, InitialState); // (phi_{i})_{x_1}
    SolutionData vForwardSpatialDerivative2(_params, _paths, _mpi, InitialState); // (phi_{i})_{x_2}
    SolutionData vForwardLaplacian(_params, _paths, _mpi, InitialState); // Lap(phi_{i})
    SolutionData vBackwardSpatialDerivative1(_params, _paths, _mpi, InitialState); // (phi^*_{i})_{x_1}
    SolutionData vBackwardSpatialDerivative2(_params, _paths, _mpi, InitialState); // (phi^*_{i})_{x_2}
    SolutionData vNonlinearTermPrevious(_params, _paths, _mpi, InitialState); // N(phi^*_{i+1})
    for (size_t i = 0; i < localGridSize; ++i) {
        vNonlinearTermPrevious[i] = complex<double>{0.0, 0.0};
    }

    for (size_t i = totalSteps; i > 0; --i) {
        dTimePoint = i * dt;
        loadForwardState(dTimePoint, i, fullSteps, stepsPerFile, vHistoryIntermediate, vHistoryRemainder, vForwardStateCurrent);

        // Laplacian(phi), phi_{x_1} and phi_{x_2} in physical space
        for (size_t j = 0; j < localGridSize; ++j) {
            vForwardSpatialDerivative1[j] = differentialOperator1[j] * vForwardStateCurrent[j];
            vForwardSpatialDerivative2[j] = differentialOperator2[j] * vForwardStateCurrent[j];
            vForwardLaplacian[j] = laplaceOperator[j] * vForwardStateCurrent[j];
        }
        _fftwPlan.ifft2InPlace(_params,vForwardSpatialDerivative1.getDataPointer());
        _fftwPlan.ifft2InPlace(_params,vForwardSpatialDerivative2.getDataPointer());
        _fftwPlan.ifft2InPlace(_params,vForwardLaplacian.getDataPointer());

        for (size_t j = 0; j < 4; ++j) {
            
            // phi^*_{x_1} and phi^*_{x_2} in Fourier space
            for (size_t k = 0; k < localGridSize; ++k) {
                vBackwardSpatialDerivative1[k] = differentialOperator1[k] * vStateCurrent[k];
                vBackwardSpatialDerivative2[k] = differentialOperator2[k] * vStateCurrent[k];
            }
            
            // phi^*, phi^*_{x_1} and phi^*_{x_2} in physical space
            _fftwPlan.ifft2InPlace(_params, vStateCurrent.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vBackwardSpatialDerivative1.getDataPointer());
            _fftwPlan.ifft2InPlace(_params, vBackwardSpatialDerivative2.getDataPointer());

            // build nonlinear term in physical space
            for (size_t k = 0; k < localGridSize; ++k) {
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
            for (size_t k = 0; k < localGridSize; ++k) {
                const complex<double> Lin = linearOperator[k];
                vNonlinearTermCurrent[k] *= dealiasingOperator[k];
                vStateNext[k] = ( (1.0 - dt * coeffBetaI[j] * Lin) * vStateCurrent[k]
                                 - dt * coeffAlphaE[j] * vNonlinearTermCurrent[k]
                                 - dt * coeffBetaE[j]  * vNonlinearTermPrevious[k] ) 
                                 / (1.0 + dt * coeffAlphaI[j] * Lin);
            }

            vStateCurrent.swapDataFrom(vStateNext);
            vNonlinearTermPrevious.swapDataFrom(vNonlinearTermCurrent);
        }
        
        // mean-zero correction
        if (localGridStart2 == 0 && localGridSize2 > 0) {
            vStateCurrent[0] = complex<double>{0.0, 0.0};
        }
    }
    vObjectiveGradient.moveDataFrom(vStateCurrent);
    setSolutionState(SolveBackwardInitialState, vObjectiveGradient);
}

void Solver::solveRiemmanianOptimization(double& dObjectiveValue, SolutionData& vObjectiveGradient, SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, 
                                     SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    
    // Finds the maximum L2 energy via Riemmanian Conjugate Gradient with Polak-Ribiere Momentum
    if (_params.getOptimizeSolution()) {
        
        // Clear old data
        if (_mpi.isRoot()) {
            filesystem::remove(_paths.getOptimizationDiagnosticsFile());
            filesystem::remove(_paths.getOptimizationLineSearchFile());
        }

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

        if (_mpi.isRoot()) {
            cout << left << setw(5)  << "Iter" << setw(18) << "Objective (J)" << setw(10) << "Adjoint" << setw(10) 
                 << "Stepsize" << setw(6) << "Brent" << setw(10) << "Forward" << setw(10) << "Objective Delta" << '\n' << string(75, '-') << '\n';
        }

        size_t iMomentumCounter = 0;
        SolutionData RawUpdate(_params, _paths, _mpi, InitialState); 
        SolutionData DirectionCurr(_params, _paths, _mpi, InitialState); 
        SolutionData DirectionPrev(_params, _paths, _mpi, InitialState); 
        SolutionData VectorTransport(_params, _paths, _mpi, InitialState);
        SolutionData ProjectedObjGradPrev(_params, _paths, _mpi, InitialState); 
        SolutionData ProjectedObjGradCurr(_params, _paths, _mpi, InitialState); 
        SolutionData TransportedProjectedObjGradPrev(_params, _paths, _mpi, InitialState); 
        SolutionData DeltaProjectedObjGrad(_params, _paths, _mpi, InitialState); 
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
        while ( abs(dObjectiveDelta) > _params.getOptimizationTolerance() && iter <= iMaxIter ) {
            if (_mpi.isRoot()) {
                cout << left << setw(5)  << iter << setw(18) << setprecision(12) << dObjectiveValue << flush;
                cout << defaultfloat << setprecision(6);
            }

            // Solve adjoint equation
            setSolutionInTime(SolveBackwardInTime, vObjectiveGradient, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            _timer.printIterationInterval(); 
            
            // Update initial data via RCG
            dObjectiveGradientSize = vObjectiveGradient.getEnergyL2(); 
            dAngleFwdICAndGradJ = vTargetStart.getInnerProductL2With(vObjectiveGradient);
            ProjectedObjGradCurr = vObjectiveGradient - (dAngleFwdICAndGradJ / _params.getInitialEnergy()) * vTargetStart;
            
            if (iter > 1) {
                dRawUpdateSize = RawUpdate.getEnergyL2();
                dAngleRawUpdateAndDirectionPrev = RawUpdate.getInnerProductL2With(DirectionPrev);                       
                dAngleRawUpdateAndProjectedObjGradPrev = RawUpdate.getInnerProductL2With(ProjectedObjGradPrev);          
                VectorTransport = (DirectionPrev - (dAngleRawUpdateAndDirectionPrev / dRawUpdateSize) * (RawUpdate)) 
                                  * sqrt(_params.getInitialEnergy() / dRawUpdateSize);          
                TransportedProjectedObjGradPrev = (ProjectedObjGradPrev - (dAngleRawUpdateAndProjectedObjGradPrev / dRawUpdateSize) * (RawUpdate))
                                  * sqrt(dManifoldSize / dRawUpdateSize); 
                DeltaProjectedObjGrad = ProjectedObjGradCurr - TransportedProjectedObjGradPrev;                       
                dAngleProjGradJAndDeltaProjGradJ = ProjectedObjGradCurr.getInnerProductL2With(DeltaProjectedObjGrad);
                dProjectedObjGradPrevSize = ProjectedObjGradPrev.getEnergyL2();     
                dMomentumSize = dAngleProjGradJAndDeltaProjGradJ / dProjectedObjGradPrevSize;                    
                iMomentumCounter = iMomentumCounter + 1;
                
                if (iter % 20 == 0) {
                    dMomentumSize = 0.0;
                }
                DirectionCurr = ProjectedObjGradCurr + (dMomentumSize * VectorTransport); 
            }
            else {
                dStepSize = 1e5;
                DirectionCurr = ProjectedObjGradCurr;
            }
                         
            dDirectionSize = DirectionCurr.getNormL2();                       
            dProjectedObjGradCurrSize = ProjectedObjGradCurr.getNormL2();
            dAscentSize = DirectionCurr.getInnerProductL2With(ProjectedObjGradCurr) / ( dDirectionSize * dProjectedObjGradCurrSize);
            if (dAscentSize < 0) {
                DirectionCurr = ProjectedObjGradCurr;     
            }        

            solveLineSearchOptimization(dStepSize, DirectionCurr, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            if (dStepSize == 0) {
                dStepSize = dAngleFwdICAndGradJ / dObjectiveGradientSize; 
            }
                                     
            RawUpdate = vTargetStart + ( dStepSize * DirectionCurr );                                   
            dRetractionSize =  sqrt(_params.getInitialEnergy()) / RawUpdate.getNormL2();                                                
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

            if (_mpi.isRoot()) {
                cout << setw(10) << dObjectiveDelta << flush << '\n';
            }        
            
            vDiagnostics = {dObjectiveValue, dObjectiveDelta, dStepSize, dManifoldSize, _timer.elapsedSeconds(), dObjectiveGradientSize, dMomentumSize};
            saveOptimizationDiagnostics(vDiagnostics);
        }

        _params.setOptimizeSolution(false);

        vHistoryIntermediate.deleteData();
        setSolutionInTime(SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
        EnergyData maxEnergyResult = getMaxEnergyL2InTimeWindow();
        vOptimalEnergySolution = { _params.getInitialEnergy(), _params.getDomainFactor1(), _params.getDomainFactor2(), _params.getTimeWindow(), dObjectiveValue, 
                                   maxEnergyResult.dTimepoint, maxEnergyResult.dEnergy };
        saveSolutionBranchAndPowerLaws(vOptimalEnergySolution);

        if (_mpi.isRoot()) {
            cout << string(75, '-') << '\n';
        }
        
        _timer.printInterval("Energy maximization problem solved at ");
    }
}

void Solver::solveLineSearchOptimization(double& dStepSize, SolutionData& DirectionCurr, SolutionData& vTargetStart, 
                                         SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    
    // Finds the optimal step size using bracketing followed by Brent's method.
    _params.setActiveLineSearch(true);

    constexpr double relativeTolerance    = 1e-5;
    constexpr double goldenExpansionFactor = 1.618034;
    constexpr double maximumExpansionFactor = 100.0;
    constexpr double goldenSectionFactor  = 0.381966;
    constexpr double numericalEpsilon   = 1e-10;
    const double initialEnergy = _params.getInitialEnergy();

    size_t iteration = 0;
    constexpr size_t maximumIterations = 1000;

    vector<double> vLineSearchHistory;
    vLineSearchHistory.reserve(maximumIterations);

    SolutionData RawUpdate(_params, _paths, _mpi, InitialState);
    SolutionData RetractedState(_params, _paths, _mpi, InitialState);
    
    // Evaluate f(step) = -J(step), where J is the terminal L2 energy. Therefore, minimizing f maximizes J.
    auto evaluateStep = [&](double step) {
        RawUpdate = vTargetStart + step * DirectionCurr;
        const double rawNorm = RawUpdate.getNormL2();

        if (rawNorm == 0.0 || !isfinite(rawNorm)) {
            return numeric_limits<double>::infinity();
        }

        const double retraction = sqrt(initialEnergy) / rawNorm;
        RetractedState = retraction * RawUpdate;
        setSolutionInTime(SolveForwardInTime, RetractedState, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
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
        double trialStepSize = middleStepSize - ((middleStepSize - nextStepSize) * interpolationTermNext 
                               - (middleStepSize - previousStepSize) * interpolationTermPrevious) / interpolationDenominator;
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
        
        if (_mpi.isRoot()) {
            cout << setw(6) << iteration << flush;
        }

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

        if (abs(bestStep - bracketMidpoint) <= doubledTolerance - 0.5 * (bracketUpperBound - bracketLowerBound)) {
            break;
        }

        bool useParabolicStep = false;
        if (abs(previousStepOffset) > absoluteTolerance) {
            const double interpolationTermPrevious = (bestStep - previousBestStep) * (bestObjective - secondPreviousBestObjective);
            double interpolationTermSecondPrevious = (bestStep - secondPreviousBestStep) * (bestObjective - previousBestObjective);
            double parabolicNumerator = (bestStep - secondPreviousBestStep) * interpolationTermSecondPrevious - (bestStep - previousBestStep) * interpolationTermPrevious;
            double parabolicDenominator = 2.0 * (interpolationTermSecondPrevious - interpolationTermPrevious);

            if (parabolicDenominator > 0.0) {
                parabolicNumerator = -parabolicNumerator;
            }

            parabolicDenominator = abs(parabolicDenominator);
            const double savedPreviousStepOffset = previousStepOffset;
            previousStepOffset = proposedStepOffset;

            if (parabolicDenominator > 0.0 && abs(parabolicNumerator) < abs(0.5 * parabolicDenominator * savedPreviousStepOffset) 
                                           && parabolicNumerator > parabolicDenominator * (bracketLowerBound - bestStep) 
                                           && parabolicNumerator < parabolicDenominator * (bracketUpperBound - bestStep)) {
                proposedStepOffset = parabolicNumerator / parabolicDenominator;
                useParabolicStep = true;

                const double proposedTrialStep = bestStep + proposedStepOffset;

                if (proposedTrialStep - bracketLowerBound < doubledTolerance || bracketUpperBound - proposedTrialStep < doubledTolerance) {
                    proposedStepOffset = copysign(absoluteTolerance, bracketMidpoint - bestStep);
                }
            }
        }

        if (!useParabolicStep) {
            previousStepOffset = (bestStep >= bracketMidpoint) ? bracketLowerBound - bestStep : bracketUpperBound - bestStep;
            proposedStepOffset = goldenSectionFactor * previousStepOffset;
        }

        const double trialStep = abs(proposedStepOffset) >= absoluteTolerance ? bestStep + proposedStepOffset : bestStep + copysign(absoluteTolerance, proposedStepOffset);
        const double trialObjective = evaluateStep(trialStep);

        if (trialObjective <= bestObjective) {
            if (trialStep >= bestStep) {
                bracketLowerBound = bestStep;
            }
            else {
                bracketUpperBound = bestStep;
            }

            secondPreviousBestStep = previousBestStep;
            secondPreviousBestObjective = previousBestObjective;
            previousBestStep = bestStep;
            previousBestObjective = bestObjective;
            bestStep = trialStep;
            bestObjective = trialObjective;
        }
        else {
            if (trialStep < bestStep) {
                bracketLowerBound = trialStep;
            }
            else {
                bracketUpperBound = trialStep;
            }

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
    _params.setActiveLineSearch(false);

    _timer.printIterationInterval();
    if (_mpi.isRoot()) {
        cout << setw(6) << iteration << flush;
    }
}

// Public functions
Solver::Solver(Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, const MPIContext& mpi) 
                : _params(params), _paths(paths), _fftwPlan(fftwPlan), _timer(timer), _mpi(mpi) {};

void Solver::setSolutionState(StateSolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        
        case SolveInitialState: {               
            if (!_params.getNumericalContinuation()) {
                setInitialCondition(vTargetState);
                
                if (_mpi.isRoot()) {
                    cout << "Generated initial guess.\n" << flush;
                }
            }
            else if (filesystem::exists(_paths.getInitialDataFile()) && _params.getOptimizeSolution()) {
                vTargetState.loadData(InitialState);
                
                if (_mpi.isRoot()) {
                    cout << "Loaded optimal initial data.\n" << flush;
                }
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
            if (_mpi.isRoot()) {
                cout << "\n";
            }    
            solveRiemmanianOptimization(dTargetValue, vObjectiveGradient, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            if (_mpi.isRoot()) {
                cout << "\n" << setprecision(12) << "Maximum Energy Amplification: " << dTargetValue << defaultfloat << setprecision(6);
            }
            break;

        case OptimizeLineSearchStepSize:
            solveLineSearchOptimization(dTargetValue, vObjectiveGradient, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            break;

    } 
    return dTargetValue;
}
