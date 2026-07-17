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

void Solver::saveSolutionBranch(const array<double, 7>& vOptimalEnergySolution) {
    vector<array<double, 7>> vSolutionsInBranch;
    
    {
        ifstream inputFile(_paths.fOptimalSolutionBranches);
        array<double, 7> solution;
        while (inputFile >> solution[0] >> solution[1] >> solution[2] >> solution[3] >> solution[4] >> solution[5]) {
            vSolutionsInBranch.push_back(solution);
        }
    }

    vSolutionsInBranch.push_back(vOptimalEnergySolution);
    sort(vSolutionsInBranch.begin(), vSolutionsInBranch.end(), [](const auto& a, const auto& b) { return a[2] < b[2]; } );
    
    ofstream outputFile(_paths.fOptimalSolutionBranches);
    outputFile << setprecision(16) << scientific;
    for (const auto& solution : vSolutionsInBranch) {
        for (size_t column = 0; column < solution.size(); ++column) {
            if (column > 0)
                outputFile << ' ';
            outputFile << solution[column];
        }
        outputFile << '\n';
    }
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

MaxEnergy Solver::getMaxEnergyInTimeWindow() {
    ifstream inputFile(_paths.fEnergyEvolution);
    double dTimepoint;
    double dEnergy;
    MaxEnergy result{ -numeric_limits<double>::infinity(), 0.0 };
    while (inputFile >> dTimepoint >> dEnergy) {
        if (dEnergy > result.dEnergy) {
            result.dEnergy = dEnergy;
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

void Solver::findContinuationForInitialData(SolutionData& vTargetState) {
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

        MaxEnergy maxEnergyResult = getMaxEnergyInTimeWindow();
        vOptimalEnergySolution = { _params.dInitialEnergy, _params.dDomainFactor1, _params.dDomainFactor2, _params.dTimeWindow, dObjectiveValue, 
                                   maxEnergyResult.dTimepoint, maxEnergyResult.dEnergy };
        saveSolutionBranch(vOptimalEnergySolution);
        cout << string(75, '-') << '\n';
        _timer.printInterval("Energy maximization problem solved at ");
    }
}

void Solver::solveLineSearchOptimization(double& dStepSize, SolutionData& DirectionCurr, SolutionData& vTargetStart, 
                                         SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    
    // Finds the optimal step size using bracketing followed by Brent's method.
    constexpr double TOL    = 1e-5;
    constexpr double GOLD   = 1.618034;
    constexpr double GLIMIT = 100.0;
    constexpr double CGOLD  = 0.381966;
    constexpr double ZEPS   = 1e-10;

    size_t iter = 0;
    constexpr size_t iMaxIter = 1000;

    vector<double> vLineSearchHistory;
    vLineSearchHistory.reserve(iMaxIter);

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
        if (iter < iMaxIter) {
            vLineSearchHistory.push_back(objectiveValue);
            ++iter;
        }
        return -objectiveValue;
    };

    // Bracketing stage
    double tA = 0.0;
    double tB = dStepSize;

    double FA = evaluateStep(tA);
    double FB = evaluateStep(tB);

    if (FB > FA) {
        swap(tA, tB);
        swap(FA, FB);
    }

    double tC = tB + GOLD * (tB - tA);
    double FC = evaluateStep(tC);

    while (FB >= FC && iter < iMaxIter) {
        const double R = (tB - tA) * (FB - FC);
        const double Q = (tB - tC) * (FB - FA);
        const double denominator = 2.0 * copysign(max(abs(Q - R), ZEPS), Q - R);
        double tP = tB - ((tB - tC) * Q - (tB - tA) * R) / denominator;
        const double tPMax = tB + GLIMIT * (tC - tB);
        double FP = 0.0;
        if ((tB - tP) * (tP - tC) > 0.0) {
            FP = evaluateStep(tP);
            if (FP < FC) {
                tA = tB;
                FA = FB;
                tB = tP;
                FB = FP;
                break;
            }
            if (FP > FB) {
                tC = tP;
                FC = FP;
                break;
            }
            tP = tC + GOLD * (tC - tB);
            FP = evaluateStep(tP);
        }
        else if ((tC - tP) * (tP - tPMax) > 0.0) {
            FP = evaluateStep(tP);
            if (FP < FC) {
                tB = tC;
                FB = FC;
                tC = tP;
                FC = FP;
                tP = tC + GOLD * (tC - tB);
                FP = evaluateStep(tP);
            }
        }
        else if ((tP - tPMax) * (tPMax - tC) >= 0.0) {
            tP = tPMax;
            FP = evaluateStep(tP);
        }
        else {
            tP = tC + GOLD * (tC - tB);
            FP = evaluateStep(tP);
        }

        tA = tB;
        tB = tC;
        tC = tP;
        FA = FB;
        FB = FC;
        FC = FP;
    }

    if (iter >= iMaxIter) {
        saveLineSearch(vLineSearchHistory);
        _timer.printIterationInterval();
        cout << setw(6) << iter << flush;
        return;
    }

    // Brent minimization stage
    double A = min(tA, tC);
    double B = max(tA, tC);
    double V = tB;
    double W = tB;
    double X = tB;
    double FX = FB;
    double FW = FB;
    double FV = FB;
    double D = 0.0;
    double E = 0.0;

    while (iter < iMaxIter) {
        const double XM   = 0.5 * (A + B);
        const double TOL1 = TOL * abs(X) + ZEPS;
        const double TOL2 = 2.0 * TOL1;
        if (abs(X - XM) <= TOL2 - 0.5 * (B - A))
            break;

        bool useParabolicStep = false;
        if (abs(E) > TOL1) {
            const double R = (X - W) * (FX - FV);
            double Q       = (X - V) * (FX - FW);
            double P       = (X - V) * Q - (X - W) * R;
            Q = 2.0 * (Q - R);
            if (Q > 0.0)
                P = -P;
            Q = abs(Q);
            const double previousE = E;
            E = D;
            if (Q > 0.0 &&
                abs(P) < abs(0.5 * Q * previousE) &&
                P > Q * (A - X) &&
                P < Q * (B - X)) {
                D = P / Q;
                useParabolicStep = true;
                const double trialPoint = X + D;
                if (trialPoint - A < TOL2 || B - trialPoint < TOL2)
                    D = copysign(TOL1, XM - X);
            }
        }

        if (!useParabolicStep) {
            E = (X >= XM) ? A - X : B - X;
            D = CGOLD * E;
        }

        const double U = abs(D) >= TOL1 ? X + D : X + copysign(TOL1, D);
        const double FU = evaluateStep(U);
        if (FU <= FX) {
            if (U >= X)
                A = X;
            else
                B = X;
            V  = W;
            FV = FW;
            W  = X;
            FW = FX;
            X  = U;
            FX = FU;
        }
        else {
            if (U < X)
                A = U;
            else
                B = U;
            if (FU <= FW || W == X) {
                V  = W;
                FV = FW;
                W  = U;
                FW = FU;
            }
            else if (FU <= FV || V == X || V == W) {
                V  = U;
                FV = FU;
            }
        }
    }
 
    dStepSize = isfinite(X) ? X : 0.0;
    saveLineSearch(vLineSearchHistory);
    _timer.printIterationInterval();
    cout << setw(6) << iter << flush;
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
