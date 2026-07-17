#ifndef PARAMETERS_H  
#define PARAMETERS_H     

#include <cmath> 
#include <string>

#include <vector>
#include <complex>
#include <array>

using namespace std;

constexpr double dPI = 3.141592653589793238462643383279502884; // reliable proxy
const complex<double> Imaginary(0.0, 1.0);

class Parameters {
private:
    vector<double> _vGridpoints1;
    vector<double> _vGridpoints2;

    vector<double> _vWavenumbersLinear1;
    vector<double> _vWavenumbersLinear2;

    vector<double> _vSpectrumNonlinear1;
    vector<double> _vSpectrumNonlinear2;
    vector<double> _vSpectrumLinear1;
    vector<double> _vSpectrumLinear2;

    vector<double> setPhysicalSpace1();
    vector<double> setPhysicalSpace2();
    vector<double> setFourierModesNonlinear1();
    vector<double> setFourierModesNonlinear2();
    vector<double> setFourierModesLinear1();
    vector<double> setFourierModesLinear2();

public:
    Parameters(int argc, char* argv[]);

    // numerical parameters
    string strInitialGuessName;  
    size_t iGridSize1;
    size_t iGridSize2;
    size_t iTotalGridSize;
    double dTimeStep;
    double dSpaceStep;
    double dInitialEnergy;
    double dEnergyFactor;
    double dDomainFactor1;
    double dDomainFactor2;
    double dDomainSize1;
    double dDomainSize2;
    double dTimeWindow;

    vector<double> vGrid1;
    vector<double> vGrid2;
    vector<double> vWavenumbersNonlinear1;
    vector<double> vWavenumbersNonlinear2;

    vector<complex<double>> vLinearOperator;
    vector<complex<double>> vLaplaceOperator;
    vector<complex<double>> vDifferentialOperator1;
    vector<complex<double>> vDifferentialOperator2;
    vector<complex<double>> vDealiasingOperator;

    vector<double> vH1Weight;
    vector<double> vH2Weight;

    vector<size_t> vRadialBin;
    size_t iMaxRadialBin;

    array<double, 4> coeffAlphaI;
    array<double, 4> coeffBetaI;
    array<double, 4> coeffAlphaE;
    array<double, 4> coeffBetaE;

    // optimization assumptions
    bool bOptimizeSolution;
    double dOptimizationTolerance;

    // initial assumptions
    bool bNumericalContinuation; 
    
    // for extended-time simulations
    double dOptimalTimeWindow; 

    size_t getIndex(size_t i, size_t j, size_t N);
    size_t iGetNumericalSteps() const;
    size_t iGetNumericalStepsPerFile() const;
    void getMathematicalOperators();

};

#endif  