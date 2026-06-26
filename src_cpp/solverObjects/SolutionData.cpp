#include "SolutionData.h"
#include "Parameters.h"

#include <complex>
#include <vector>

#include <iostream>  

using Complex = std::complex<double>;

SolutionData::SolutionData(const Parameters& params, SolutionDataType storedDataType) {
    int iTotalGridPoints = params.iGridSize1 * params.iGridSize2;

    switch (storedDataType) {
        case InitialState:
        case TerminalState:
        case BackwardInitialState:
            vData.resize(iTotalGridPoints * 1);
            std::cout << "Gridpoints: " << vData.size() << "\n";
            break;
        case IntermediateHistory:   
            vData.resize(iTotalGridPoints * params.iGetNumericalStepsPerFile() );
            std::cout << "Intermediate States: " << vData.size() / iTotalGridPoints << "\n";
            break;
        case RemainderHistory:
            vData.resize(iTotalGridPoints * ( params.iGetNumericalSteps() % params.iGetNumericalStepsPerFile()) );
            std::cout << "Remainder States: " << vData.size() / iTotalGridPoints << "\n";
            break;
    }
}

std::size_t SolutionData::getSize() const {
    return vData.size();
}

std::vector<Complex>& SolutionData::getData() {
    return vData;
}

Complex* SolutionData::getDataPointer() {
    return vData.data();
}

/*
const Complex *SolutionData::getDataPointer() const
{
    return vec_Data.data();
}
*/

/*
std::vector<Complex> &SolutionData::vec_SetStateInHistory(const std::vector<Complex> &state)
{
    vec_StateHistory = state;
}

void &SolutionData::saveHistory(const std::vector<Complex> &state)
{
    vec_StateHistory = state;
}

void &SolutionData::saveInitialState(const std::vector<Complex> &state)
{
    vec_StateHistory = state;
}

void &SolutionData::saveTerminalState(const std::vector<Complex> &state)
{
    vec_StateHistory = state;
}

*/