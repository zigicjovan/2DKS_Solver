#ifndef PARAMETERS_H  
#define PARAMETERS_H     

#include <cmath> 
#include <string>

inline constexpr double dPI = 3.141592653589793238462643383279502884; // reliable proxy

class Parameters {
private:

public:
    // numerical parameters
    std::string strInitalGuessName;  
    int iGridSize1;
    int iGridSize2;
    double dTimeStep;
    double dInitialEnergy;
    double dDomainFactor1;
    double dDomainFactor2;
    double dDomainSize1;
    double dDomainSize2;
    double dTimeWindow;

    // optimization assumptions
    bool bOptimizeSolution;
    double dOptimizationTolerance;

    // initial assumptions
    bool bNumericalContinuation; 
    
    // for extended-time simulations
    double dOptimalTimeWindow; 

    int iGetNumericalSteps() const {   
        return static_cast<int>(std::round(dTimeWindow / dTimeStep));
    }

    // fixed file size
    int iGetNumericalStepsPerFile() const {   
        return static_cast<int>(std::round( 4e7 / (iGridSize1 * iGridSize2) ));
    }
};

#endif  