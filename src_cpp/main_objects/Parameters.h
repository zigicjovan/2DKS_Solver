// what: declare prototype for chosen numerical parameters and testing assumptions

#ifndef PARAMETERS_H  
#define PARAMETERS_H     

#include <cmath> // why: for std::round
#include <string>

inline constexpr double d_PI = 3.141592653589793238462643383279502884; // why: reliable proxy for pi value

struct Parameters
{
    // numerical parameters
    std::string str_InitalGuessName;
    
    int n_GridSize_1;
    int n_GridSize_2;
    double d_TimeStep;

    double d_InitialEnergy;
    double d_DomainFactor_1;
    double d_DomainFactor_2;
    double d_DomainSize_1;
    double d_DomainSize_2;
    double d_TimeWindow;

    int n_numericalSteps() const 
    {   // how: divide time window by time step to approximate number of required steps
        return static_cast<int>(std::round(d_TimeWindow / d_TimeStep));
    }

    int n_numericalStepsPerFile() const 
    {   // why: maintain fixed directory forward data file size
        return static_cast<int>(std::round( 4e7 / (n_GridSize_1 * n_GridSize_2) ));
    }

    // optimization assumptions
    bool b_OptimizeSolution;
    double d_OptimizationTolerance;

    // initial assumptions
    bool b_NumericalContinuation; // why: use nearby solution as an approximate initial guess
    double d_OptimalTimeWindow; // why: for extended-time simulations

};

#endif  