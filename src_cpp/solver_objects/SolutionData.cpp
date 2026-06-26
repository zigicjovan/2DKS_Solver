#include "SolutionData.h"
#include "Parameters.h"

#include <complex>
#include <vector>

#include <iostream>  

using Complex = std::complex<double>;

SolutionData::SolutionData(const Parameters &params, SolutionDataType storedDataType)
{
    int n_TotalGridPoints = params.n_GridSize_1 * params.n_GridSize_2;

    switch (storedDataType)
    {
        case SolutionDataType::InitialState:
        case SolutionDataType::TerminalState:
        case SolutionDataType::BackwardInitialState:
        {   
            vec_Data.resize(n_TotalGridPoints * 1);
            std::cout   << "Gridpoints: "           << vec_Data.size() << "\n";
            break;
        }
        case SolutionDataType::IntermediateHistory:   
        {
            vec_Data.resize(n_TotalGridPoints * params.n_numericalStepsPerFile() );
            std::cout   << "Intermediate States: "  << vec_Data.size() / n_TotalGridPoints << "\n";
            break;
        }
        case SolutionDataType::RemainderHistory:
        {   
            vec_Data.resize(n_TotalGridPoints * ( params.n_numericalSteps() % params.n_numericalStepsPerFile()) );
            std::cout   << "Remainder States: "     << vec_Data.size() / n_TotalGridPoints << "\n";
            break;
        }
    }
}

std::size_t SolutionData::getSize() const 
{
    return vec_Data.size();
}

std::vector<Complex> &SolutionData::getData()
{
    return vec_Data;
}

Complex *SolutionData::getDataPointer()
{
    return vec_Data.data();
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