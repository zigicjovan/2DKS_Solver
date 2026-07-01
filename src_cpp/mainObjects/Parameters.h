#ifndef PARAMETERS_H  
#define PARAMETERS_H     

#include <cmath> 
#include <string>

#include <vector>
#include <complex>
#include <array>

using Complex = std::complex<double>;
constexpr double dPI = 3.141592653589793238462643383279502884; // reliable proxy
const std::complex<double> Imaginary(0.0, 1.0);
std::size_t getIndex(std::size_t i, std::size_t j, std::size_t N);

class Parameters {
private:
    std::vector<double> _vGridpoints1;
    std::vector<double> _vGridpoints2;

    std::vector<double> _vWavenumbersLinear1;
    std::vector<double> _vWavenumbersLinear2;

    std::vector<double> _vSpectrumNonlinear1;
    std::vector<double> _vSpectrumNonlinear2;
    std::vector<double> _vSpectrumLinear1;
    std::vector<double> _vSpectrumLinear2;

    std::vector<double> setPhysicalSpace1();
    std::vector<double> setPhysicalSpace2();
    std::vector<double> setFourierModesNonlinear1();
    std::vector<double> setFourierModesNonlinear2();
    std::vector<double> setFourierModesLinear1();
    std::vector<double> setFourierModesLinear2();

public:
    // numerical parameters
    std::string strInitialGuessName;  
    std::size_t iGridSize1;
    std::size_t iGridSize2;
    std::size_t iTotalGridSize;
    double dTimeStep;
    double dSpaceStep;
    double dInitialEnergy;
    double dDomainFactor1;
    double dDomainFactor2;
    double dDomainSize1;
    double dDomainSize2;
    double dTimeWindow;

    std::vector<double> vGrid1;
    std::vector<double> vGrid2;
    std::vector<double> vWavenumbersNonlinear1;
    std::vector<double> vWavenumbersNonlinear2;

    std::vector<Complex> vLinearOperator;
    std::vector<Complex> vLaplaceOperator;
    std::vector<Complex> vDifferentialOperator1;
    std::vector<Complex> vDifferentialOperator2;

    std::array<double, 4> coeffAlphaI;
    std::array<double, 4> coeffBetaI;
    std::array<double, 4> coeffAlphaE;
    std::array<double, 4> coeffBetaE;

    // optimization assumptions
    bool bOptimizeSolution;
    double dOptimizationTolerance;

    // initial assumptions
    bool bNumericalContinuation; 
    
    // for extended-time simulations
    double dOptimalTimeWindow; 

    int iGetNumericalSteps() const;
    int iGetNumericalStepsPerFile() const;
    void getMathematicalOperators();

};

#endif  