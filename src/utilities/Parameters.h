#ifndef PARAMETERS_H  
#define PARAMETERS_H     

#include "MPIContext.h"

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
    
    // grid parameters
    vector<double> _vGridpoints1;
    vector<double> _vGridpoints2;

    vector<double> _vWavenumbersLinear1;
    vector<double> _vWavenumbersLinear2;

    vector<double> _vSpectrumNonlinear1;
    vector<double> _vSpectrumNonlinear2;
    vector<double> _vSpectrumLinear1;
    vector<double> _vSpectrumLinear2;

    vector<double> _vGrid1;
    vector<double> _vGrid2;
    vector<double> _vWavenumbersNonlinear1;
    vector<double> _vWavenumbersNonlinear2;

    // solver parameters
    string _strInitialGuessName;  
    size_t _iGridSize1;
    size_t _iGridSize2;
    size_t _iTotalGridSize;
    size_t _iSavedStates;

    double _dTimeStep;
    double _dSpaceStep;
    double _dInitialEnergy;
    double _dEnergyFactor;
    double _dDomainFactor1;
    double _dDomainFactor2;
    double _dDomainSize1;
    double _dDomainSize2;
    double _dTimeWindow;

    vector<complex<double>> _vLinearOperator;
    vector<complex<double>> _vLaplaceOperator;
    vector<complex<double>> _vDifferentialOperator1;
    vector<complex<double>> _vDifferentialOperator2;
    vector<unsigned char> _vDealiasingOperator;

    vector<double> _vH1Weight;
    vector<double> _vH2Weight;
    vector<size_t> _vRadialBin;
    size_t _iMaxRadialBin;
    double _dRadialBinWidth;

    array<double, 4> _coeffAlphaI;
    array<double, 4> _coeffBetaI;
    array<double, 4> _coeffAlphaE;
    array<double, 4> _coeffBetaE;

    // optimization assumptions
    bool _bOptimizeSolution;
    double _dOptimizationTolerance;
    bool _bActiveLineSearch;

    // initial assumptions
    bool _bNumericalContinuation; 
    
    // for extended-time simulations
    double _dOptimalTimeWindow; 

    vector<double> setPhysicalSpace1();
    vector<double> setPhysicalSpace2();
    vector<double> setFourierModesNonlinear1();
    vector<double> setFourierModesNonlinear2();
    vector<double> setFourierModesLinear1();
    vector<double> setFourierModesLinear2();

public:
    Parameters(int argc, char* argv[]);

    void getMathematicalOperators(const MPIContext& mpi);

    size_t getIndex(size_t i, size_t j, size_t N);
    size_t getNumericalSteps() const;
    size_t getNumericalStepsPerFile() const;

    const vector<double>& getGrid1() const;
    const vector<double>& getGrid2() const;
    const vector<double>& getWavenumbersNonlinear1() const;
    const vector<double>& getWavenumbersNonlinear2() const;

    const string& getInitialGuessName() const;
    size_t getGridSize1() const;
    size_t getGridSize2() const;
    size_t getTotalGridSize() const;
    size_t getSavedStates() const;

    double getTimeStep() const;
    double getSpaceStep() const;
    double getInitialEnergy() const;
    double getEnergyFactor() const;
    double getDomainFactor1() const;
    double getDomainFactor2() const;
    double getDomainSize1() const;
    double getDomainSize2() const;
    double getTimeWindow() const;

    const vector<complex<double>>& getLinearOperator() const;
    const vector<complex<double>>& getLaplaceOperator() const;
    const vector<complex<double>>& getDifferentialOperator1() const;
    const vector<complex<double>>& getDifferentialOperator2() const;
    const vector<unsigned char>& getDealiasingOperator() const;

    const vector<double>& getH1Weight() const;
    const vector<double>& getH2Weight() const;
    const vector<size_t>& getRadialBin() const;
    size_t getMaxRadialBin() const;
    double getRadialBinWidth() const;

    const array<double, 4>& getCoeffAlphaI() const;
    const array<double, 4>& getCoeffBetaI() const;
    const array<double, 4>& getCoeffAlphaE() const;
    const array<double, 4>& getCoeffBetaE() const;

    bool getOptimizeSolution() const;
    void setOptimizeSolution(bool optimizeSolution);
    bool getActiveLineSearch() const;
    void setActiveLineSearch(bool optimizeSolution);
    double getOptimizationTolerance() const;
    bool getNumericalContinuation() const;
    double getOptimalTimeWindow() const;

};

#endif  