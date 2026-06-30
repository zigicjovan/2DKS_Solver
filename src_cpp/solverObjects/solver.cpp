#include "solver.h"
#include "Parameters.h"
#include "Pathnames.h"
#include "SolutionData.h"
#include "updateDirectory.h"

#include <filesystem>
#include <optional>
#include <regex>
#include <string>
#include <iostream>

std::filesystem::path continuedFile;
double bestT = -1.0;

void setSolutionState(const Parameters& params, Pathnames& paths, SolutionType targetType, SolutionData& vTargetState) {
    switch (targetType) {
        case SolveInitialState: {          
       
            if (params.bNumericalContinuation == 0) {
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

                    // Does this file belong to this testcase?
                    if (filename.find(paths.strTestcaseGeneric.str()) == std::string::npos)
                        continue;

                    std::size_t pos1 = filename.find("_T_"); // Find "_T_"
                    if (pos1 == std::string::npos)
                        continue;
                    pos1 += 3; // skip "_T_"
                    
                    std::size_t pos2 = filename.find("_opt_", pos1); // Find "_opt_"
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
                        std::cout << "No saved data, generated initial guess.\n";
                    }
            }
            break;
        }
        case SolveTerminalState: {
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break;
    }  
}

void setSolutionInTime(const Parameters& params, const Pathnames& paths, SolutionType targetType, SolutionData& vTargetStart, 
    SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case SolveForwardInTime: {
            break;
        }
        case SolveBackwardInTime: {
            break;
        }
        case SolveInitialState:
        case SolveTerminalState:
        case OptimizeEnergyAmplification:
        case OptimizeLineSearchStepSize:
            break;
    }  
}

double getOptimalSolution(const Parameters& params, const Pathnames& paths, SolutionType targetType, SolutionData& vObjectiveGradient, 
    SolutionData& vTargetStart, SolutionData& vHistoryIntermediate, SolutionData& vHistoryRemainder, SolutionData& vTargetEnd) {
    switch (targetType) {
        case OptimizeEnergyAmplification: {
            break;
        }
        case OptimizeLineSearchStepSize: {
            break;
        }
        case SolveForwardInTime:
        case SolveBackwardInTime:
        case SolveInitialState:
        case SolveTerminalState:
            break;
    } 
    return 0.0;
}
