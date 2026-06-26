#include "Parameters.h"
#include "readParameterSettings.h"

#include <iostream>  
#include <string>

Parameters readParameterSettings(int argc, char* argv[])
{
    Parameters params;

    params.str_InitalGuessName      = argv[1];

    params.n_GridSize_1             = std::stoi(argv[2]);
    params.n_GridSize_2             = std::stoi(argv[3]);
    params.d_TimeStep               = std::stod(argv[4]);

    params.d_InitialEnergy          = std::pow(10.0, std::stod(argv[5]));
    params.d_DomainFactor_1         = std::stod(argv[6]);
    params.d_DomainFactor_2         = std::stod(argv[7]);
    params.d_DomainSize_1           = 2 * d_PI * std::stod(argv[6]);
    params.d_DomainSize_2           = 2 * d_PI * std::stod(argv[7]);
    params.d_TimeWindow             = std::pow(10.0, std::stod(argv[8]));

    params.b_OptimizeSolution       = std::stoi(argv[9]);
    params.d_OptimizationTolerance  = std::stod(argv[10]);
    params.b_NumericalContinuation  = std::stoi(argv[11]);
    if (params.b_OptimizeSolution == true)
    {
        params.d_OptimalTimeWindow  = params.d_TimeWindow;
    }
    else
    {
        params.d_OptimalTimeWindow = std::pow(10.0, std::stod(argv[12]));
    }

    std::cout   << "Parameter settings:\nIC "     << params.str_InitalGuessName 
                << ", N_x1 "                      << params.n_GridSize_1 
                << ", N_x2 "                      << params.n_GridSize_2       
                << ", dt "                        << params.d_TimeStep         
                << ", K "                         << params.d_InitialEnergy    
                << ", ell1 "                      << params.d_DomainFactor_1     
                << ", ell2 "                      << params.d_DomainFactor_2  
                << ", L_1 "                       << params.d_DomainSize_1     
                << ", L_2 "                       << params.d_DomainSize_2     
                << ", T "                         << params.d_TimeWindow       
                << ", nsteps "                    << params.n_numericalSteps() 
                << ", nstepsFile "                << params.n_numericalStepsPerFile() 
                << ", opt "                       << params.b_OptimizeSolution       
                << ", tol "                       << params.d_OptimizationTolerance  
                << ", cont "                      << params.b_NumericalContinuation  
                << ", optT "                      << params.d_OptimalTimeWindow      
                << std::endl;

    return params;
}