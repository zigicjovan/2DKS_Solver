#include "solver.h"
#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "updateDirectory.h"
#include "FFTWPlanner.h"
#include "Timer.h"

#include <filesystem>
#include <string>
#include <iostream>

#include <complex>
#include <cstddef>
#include <random>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using Complex = std::complex<double>;

void setInitialCondition(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, SolutionData& vTargetState) {  
    if (params.strInitialGuessName == "s1") {
        int stateNumber = 0;
        for (std::size_t j = 0; j < params.iGridSize2; ++j) {
            for (std::size_t i = 0; i < params.iGridSize1; ++i) {
                const std::size_t p = getIndex(j, i, params.iGridSize1);
                const double x = params.vGrid1[p];
                const double y = params.vGrid2[p];
                vTargetState(params, i, j, stateNumber) = Complex(
                    std::sin(x / params.dDomainFactor1 + y / params.dDomainFactor2)
                    + std::sin(x / params.dDomainFactor1) + std::sin(y / params.dDomainFactor2), 0.0);
            }
        }
        Complex* vTempState = vTargetState.getStatePointer(params, stateNumber);
        fftwPlan.fft2InPlace(vTempState);
    }
    else if (params.strInitialGuessName == "randfour") {
        const int stateNumber = 0;

        const double kmax  = std::max(params.iGridSize1, params.iGridSize2) / 2.0;
        const double delta = 6.0 / kmax;

        std::mt19937 rng(std::random_device{}());
        std::uniform_real_distribution<double> dist(0.0, 2.0 * M_PI);

        Complex* vTempState = vTargetState.getStatePointer(params, stateNumber);

        // Build randomized Fourier coefficients
        for (std::size_t j = 0; j < params.iGridSize2; ++j) {
            for (std::size_t i = 0; i < params.iGridSize1; ++i) {
                const std::size_t p = getIndex(j, i, params.iGridSize1);
                const double k1 = params.vWavenumbersNonlinear1[i];
                const double k2 = params.vWavenumbersNonlinear2[j];
                const double kabs = std::sqrt(k1 * k1 + k2 * k2);
                const double A = std::exp(-delta * kabs);
                const double phi = dist(rng);
                vTempState[p] = A * Complex(std::cos(phi), std::sin(phi));
            }
        }
        vTempState[0] = Complex(0.0, 0.0); // zero mean mode
        fftwPlan.ifft2InPlace(params, vTempState);  // real(ifft2(uhat))
        for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
            vTempState[p] = Complex(vTempState[p].real(), 0.0);
        fftwPlan.fft2InPlace(vTempState);
    }
    else {
        throw std::runtime_error(
            "Unknown initial condition: " + params.strInitialGuessName);
    }
    vTargetState.setInitialEnergyL2(params);
}

void setSolutionState(const Parameters& params, Pathnames& paths, FFTWPlanner& fftwPlan, SolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        case SolveInitialState: {               
            if (params.bNumericalContinuation == 0) {
                setInitialCondition(params, paths, fftwPlan, vTargetState);
                std::cout << "Generated initial guess.\n";
            }
            else if (std::filesystem::exists(paths.fOptimalInitialData)) {
                loadData(paths, vTargetState, InitialState);
                std::cout << "Loaded optimal initial data.\n";
            }
            else {
                std::filesystem::path fContinuedFile;
                std::filesystem::path fTempFile;
                double dBestT = params.dTimeWindow + 1.0;            
                
                // Search for nearest previous time window in directory
                for (const auto& entry : std::filesystem::directory_iterator(paths.dirOptimalInitialData)) {
                    if (!entry.is_regular_file())
                        continue;
                    std::string filename = entry.path().filename().string();
                    if (filename.find(paths.strTestcaseGeneric.str()) == std::string::npos)
                        continue;
                    std::size_t pos1 = filename.find("_T_"); 
                    if (pos1 == std::string::npos)
                        continue;
                    pos1 += 3; // skip "_T_"
                    std::size_t pos2 = filename.find("_opt_", pos1); 
                    if (pos2 == std::string::npos)
                        continue;
                    double T = std::stod(filename.substr(pos1, pos2 - pos1)); // Extract T
                    if (T <= params.dTimeWindow && ( T > dBestT || dBestT > params.dTimeWindow)) {
                        dBestT = T;
                        fContinuedFile = entry.path();
                    }
                }
                if (dBestT <= params.dTimeWindow) {
                    fTempFile = paths.fOptimalInitialData;
                    paths.fOptimalInitialData = fContinuedFile;
                    loadData(paths, vTargetState, InitialState);
                    std::cout << "Loaded continued initial data from T = " << dBestT << "\n" ;
                    paths.fOptimalInitialData = fTempFile; 
                }
                else {
                    setInitialCondition(params, paths, fftwPlan, vTargetState);
                    std::cout << "No saved data, generated initial guess.\n";
                }
            }
            saveData(paths, vTargetState, InitialState);
            break;
        }
        case SolveTerminalState: {
            // TO DO: take final state from vHistoryRemainder
            saveData(paths, vTargetState, TerminalState);
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break; // Do nothing
    }  
}

void setSolutionInTime(const Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, SolutionType targetType, SolutionData& vTargetStart, 
    SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case SolveForwardInTime: {
            // TO DO: solve forward problem using IMEX RK4 [create function]
            // saveData: IntermediateHistory, RemainderHistory, TerminalState
            // saveData(paths, vHistoryIntermediate, IntermediateHistory, params.dTimeWindow);
            // saveData(paths, vHistoryRemainder, RemainderHistory, dCurrentT);
            if (params.bOptimizeSolution == 0) {
                // save EnergyEvolution, FourierSpectrumEvolution
            }
            timer.printInterval("");
            break;
        }
        case SolveBackwardInTime: {
            // TO DO: solve backward problem using IMEX RK4 [create function]
            // loadData(paths, vHistoryIntermediate, IntermediateHistory, params.dTimeWindow);
            // loadData(paths, vHistoryRemainder, RemainderHistory, dCurrentT);
            // saveData: BackwardInitialState
            timer.printInterval("");
            break;
        }
        case SolveInitialState:
        case SolveTerminalState:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break; // Do nothing
    }  
}

double getOptimalSolution(Parameters& params, const Pathnames& paths, FFTWPlanner& fftwPlan, Timer& timer, SolutionType targetType, SolutionData& vObjectiveGradient, 
    SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case OptimizeEnergyAmplification: {
            if (params.bOptimizeSolution == 1) {
                // TO DO: solve RCG  [create function]
                // saveData(paths, vTargetStart, InitialState);
                // save SolutionBranch, OptDiagnostics
                
                // save diagnostics after completion
                params.bOptimizeSolution = 0; 
                setSolutionInTime(params, paths, fftwPlan, timer, SolveForwardInTime, vTargetStart, vHistoryIntermediate, vHistoryRemainder, vTargetEnd);
            }
            timer.printInterval("");
            break;
        }
        case OptimizeLineSearchStepSize: {
            // TO DO: solve Brent [create function]
            // save OptLineSearch
            timer.printInterval("");
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case SolveInitialState:
        case SolveTerminalState:
            break; // Do nothing
    } 
    return 0.0;
}
