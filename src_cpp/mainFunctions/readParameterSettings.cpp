#include "Parameters.h"
#include "readParameterSettings.h"

#include <iostream>  
#include <string>

Parameters readParameterSettings(int argc, char* argv[]) {
    Parameters params;

    params.strInitialGuessName = argv[1];
    params.iGridSize1 = std::stoi(argv[2]);
    params.iGridSize2 = std::stoi(argv[3]);
    params.iTotalGridSize = params.iGridSize1 * params.iGridSize2;
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

    params.vLinearOperator.resize(params.iTotalGridSize);
    params.vLaplaceOperator.resize(params.iTotalGridSize);
    params.vDifferentialOperator1.resize(params.iTotalGridSize);
    params.vDifferentialOperator2.resize(params.iTotalGridSize);
    params.getMathematicalOperators();

    params.coeffAlphaI = { 343038331393.0 / 1130875731271.0,
                       288176579239.0 / 1140253497719.0,
                       253330171251.0 / 677500478386.0,
                       189462239225.0 / 1091147436423.0 };
    params.coeffBetaI = { 35965327958.0 / 140127563663.0,
                      19632212512.0 / 2700543775099.0,
                     -173747147147.0 / 351772688865.0,
                      91958533623.0 / 727726057489.0 };
    params.coeffAlphaE = { 14.0 / 25.0,
                       777974228744.0 / 1346157007247.0,
                       251277807242.0 / 1103637129625.0,
                       113091689455.0 / 220187950967.0 };
    params.coeffBetaE = { 0.0,
                     -251352885992.0 / 790610919619.0,
                     -383714262797.0 / 1103637129625.0,
                     -403360439203.0 / 1888264787188.0 };
    
    if (params.bOptimizeSolution == true)
        params.dOptimalTimeWindow = params.dTimeWindow;
    else
        params.dOptimalTimeWindow = std::pow(10.0, std::stod(argv[12]));

    std::cout << "Parameter settings:\nIC " << params.strInitialGuessName 
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