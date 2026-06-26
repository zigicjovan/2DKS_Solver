#include "Parameters.h"
#include "readParameterSettings.h"

#include <iostream>  
#include <string>

Parameters readParameterSettings(int argc, char* argv[]) {
    Parameters params;

    params.strInitalGuessName = argv[1];
    params.iGridSize1 = std::stoi(argv[2]);
    params.iGridSize2 = std::stoi(argv[3]);
    params.dTimeStep = std::stod(argv[4]);
    params.dInitialEnergy = std::pow(10.0, std::stod(argv[5]));
    params.dDomainFactor1 = std::stod(argv[6]);
    params.dDomainFactor2 = std::stod(argv[7]);
    params.dDomainSize1 = 2 * dPI * std::stod(argv[6]);
    params.dDomainSize2 = 2 * dPI * std::stod(argv[7]);
    params.dTimeWindow = std::pow(10.0, std::stod(argv[8]));
    params.bOptimizeSolution = std::stoi(argv[9]);
    params.dOptimizationTolerance = std::stod(argv[10]);
    params.bNumericalContinuation = std::stoi(argv[11]);
    
    if (params.bOptimizeSolution == true)
        params.dOptimalTimeWindow = params.dTimeWindow;
    else
        params.dOptimalTimeWindow = std::pow(10.0, std::stod(argv[12]));

    std::cout << "Parameter settings:\nIC " << params.strInitalGuessName 
              << ", N_x1 " << params.iGridSize1 
              << ", N_x2 " << params.iGridSize2       
              << ", dt " << params.dTimeStep         
              << ", K " << params.dInitialEnergy    
              << ", ell1 " << params.dDomainFactor1     
              << ", ell2 " << params.dDomainFactor2  
              << ", L_1 " << params.dDomainSize1     
              << ", L_2 " << params.dDomainSize2     
              << ", T " << params.dTimeWindow       
              << ", nsteps " << params.iGetNumericalSteps() 
              << ", nstepsFile " << params.iGetNumericalStepsPerFile() 
              << ", opt " << params.bOptimizeSolution       
              << ", tol " << params.dOptimizationTolerance  
              << ", cont " << params.bNumericalContinuation  
              << ", optT " << params.dOptimalTimeWindow      
              << std::endl;

    return params;
}