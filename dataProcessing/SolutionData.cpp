#include "SolutionData.h"
#include "Parameters.h"

#include <complex>
#include <vector> 
#include <filesystem>
#include <iostream> 

using namespace std;

// private functions
void SolutionData::writeDistributedFile(const filesystem::path& filename) const {
    const size_t localAllocationSize = _mpi.getLocalAllocationSize();
    const size_t localGridSize = _mpi.getLocalGridSize();
    const size_t globalGridSize = _params.getTotalGridSize();
    const size_t globalGridOffset = _mpi.getLocalGridStart2() * _params.getGridSize1();

    if (localAllocationSize == 0) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const size_t stateCount = _vData.size() / localAllocationSize;
    const size_t localByteCount = localGridSize * sizeof(complex<double>);

    if (localByteCount > static_cast<size_t>(INT_MAX)) {
        if (_mpi.isRoot()) {
            cerr << "ERROR: Local MPI-IO write exceeds the supported MPI count.\n";
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_File file;
    const int openResult = MPI_File_open(MPI_COMM_WORLD, filename.string().c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &file);

    if (openResult != MPI_SUCCESS) {
        if (_mpi.isRoot()) {
            cerr << "ERROR: Could not open MPI output file: " << filename << '\n';
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const MPI_Offset completeFileSize = static_cast<MPI_Offset>(stateCount) * static_cast<MPI_Offset>(globalGridSize) * static_cast<MPI_Offset>(sizeof(complex<double>));
    MPI_File_set_size(file, completeFileSize);

    for (size_t state = 0; state < stateCount; ++state) {
        const MPI_Offset fileOffset = ( static_cast<MPI_Offset>(state) * static_cast<MPI_Offset>(globalGridSize) + static_cast<MPI_Offset>(globalGridOffset) )
                                      * static_cast<MPI_Offset>(sizeof(complex<double>));
        MPI_File_write_at_all(file, fileOffset, const_cast<complex<double>*>(getStatePointer(state)), static_cast<int>(localByteCount), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&file);
}

void SolutionData::readDistributedFile(const filesystem::path& filename) {
    const size_t localAllocationSize = _mpi.getLocalAllocationSize();
    const size_t localGridSize = _mpi.getLocalGridSize();
    const size_t globalGridSize = _params.getTotalGridSize();
    const size_t globalGridOffset = _mpi.getLocalGridStart2() * _params.getGridSize1();

    if (localAllocationSize == 0) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    const size_t stateCount = _vData.size() / localAllocationSize;
    const size_t localByteCount = localGridSize * sizeof(complex<double>);

    if (localByteCount > static_cast<size_t>(INT_MAX)) {
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    MPI_File file;
    const int openResult = MPI_File_open(MPI_COMM_WORLD, filename.string().c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &file);

    if (openResult != MPI_SUCCESS) {
        if (_mpi.isRoot()) {
            cerr << "ERROR: Could not open MPI input file: " << filename << '\n';
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    for (size_t state = 0; state < stateCount; ++state) {
        const MPI_Offset fileOffset = ( static_cast<MPI_Offset>(state) * static_cast<MPI_Offset>(globalGridSize) + static_cast<MPI_Offset>(globalGridOffset) )
                                      * static_cast<MPI_Offset>(sizeof(complex<double>));
        MPI_File_read_at_all(file, fileOffset, getStatePointer(state), static_cast<int>(localByteCount), MPI_BYTE, MPI_STATUS_IGNORE);
    }

    MPI_File_close(&file);
}

// public functions
filesystem::path SolutionData::appendTimeStep(const filesystem::path& path, double dCurrentT) {
    return path.parent_path() / (path.stem().string() + "_" + to_string(dCurrentT) + path.extension().string());
}

SolutionData::SolutionData(const Parameters& params, const Pathnames& paths, const MPIContext& mpi, SolutionDataType storedDataType) 
                           : _params(params), _paths(paths), _mpi(mpi) {
    const size_t stateStorageSize = _mpi.getLocalAllocationSize();

    switch (storedDataType) {
        
        case InitialState: // fall through
        
        case TerminalState: // fall through
        
        case BackwardInitialState:
            _vData.resize(stateStorageSize);
            break;
    
        case IntermediateHistory:   
            _vData.resize(stateStorageSize * _params.getNumericalStepsPerFile() );
            break;
        
        case RemainderHistory:
            _vData.resize(stateStorageSize * ( _params.getNumericalSteps() % _params.getNumericalStepsPerFile()) );
            break;
    }
}

size_t SolutionData::getSize() const {
    return _vData.size();
}

void SolutionData::getData(SolutionData& stateData, size_t stateIndex) const {
    const size_t stateStorageSize = _mpi.getLocalAllocationSize();
    const size_t iOffset = stateIndex * stateStorageSize;
    copy(_vData.begin() + iOffset, _vData.begin() + iOffset + stateStorageSize, stateData._vData.begin());
}

void SolutionData::setData(const SolutionData& stateData, size_t stateIndex) {
    const size_t stateStorageSize = _mpi.getLocalAllocationSize();
    const size_t iOffset = stateIndex * stateStorageSize;
    copy(stateData._vData.begin(), stateData._vData.begin() + stateStorageSize, _vData.begin() + iOffset);
}

SolutionData& SolutionData::operator=(const SolutionData& otherData) {
    if (this != &otherData) {
        _vData = otherData._vData;
    }
    return *this;
}

SolutionData& SolutionData::operator+=(const SolutionData& otherData) {
    for (size_t i = 0; i < _vData.size(); ++i) {
        _vData[i] += otherData._vData[i];
    }
    return *this;
}

SolutionData& SolutionData::operator-=(const SolutionData& otherData) {
    for (size_t i = 0; i < _vData.size(); ++i) {
        _vData[i] -= otherData._vData[i];
    }
    return *this;
}

SolutionData& SolutionData::operator*=(double scalar) {
    for (complex<double>& value : _vData) {
        value *= scalar;
    }
    return *this;
}

SolutionData operator+(SolutionData lhs, const SolutionData& rhs) {
    lhs += rhs;
    return lhs;
}

SolutionData operator-(SolutionData lhs, const SolutionData& rhs) {
    lhs -= rhs;
    return lhs;
}

SolutionData operator*(SolutionData data, double scalar) {
    data *= scalar;
    return data;
}

SolutionData operator*(double scalar, SolutionData data) {
    data *= scalar;
    return data;
}

complex<double>& SolutionData::operator()(size_t i, size_t j, size_t stateNumber) {
    const size_t stateStorageSize = _mpi.getLocalAllocationSize();
    return _vData[stateNumber * stateStorageSize + j * _params.getGridSize1() + i];
}

complex<double>& SolutionData::operator[](size_t idx) {
    return _vData[idx];
}

complex<double>& SolutionData::atState(size_t p, size_t stateNumber) {
    return _vData[stateNumber * _mpi.getLocalAllocationSize() + p];
}

complex<double>* SolutionData::getStatePointer(size_t stateNumber) {
    return _vData.data() + stateNumber * _mpi.getLocalAllocationSize();
}

const complex<double>* SolutionData::getStatePointer(size_t stateNumber) const {
    return _vData.data() + stateNumber * _mpi.getLocalAllocationSize();
}

complex<double>* SolutionData::getDataPointer() {
    return _vData.data();
}

void SolutionData::swapDataFrom(SolutionData& otherData) noexcept {
    _vData.swap(otherData._vData);
}

void SolutionData::moveDataFrom(SolutionData& otherData) noexcept {
    _vData = move(otherData._vData);
}

void SolutionData::setInitialEnergyL2() {
    double energyL2 = getEnergyL2();
    double normL2 = sqrt(energyL2);
    double scale = sqrt(_params.getInitialEnergy()) / normL2;
    const size_t localGridSize = _mpi.getLocalGridSize();

    for (size_t i = 0; i < localGridSize; ++i) {
        _vData[i] *= scale;
    }
}

double SolutionData::getInnerProductL2With(SolutionData& otherData) const {
    double localInnerProduct = 0.0;
    const size_t localGridSize = _mpi.getLocalGridSize();

    for (size_t i = 0; i < localGridSize; ++i) {
        localInnerProduct += real(conj(_vData[i]) * otherData._vData[i]);
    }

    double globalInnerProduct = 0.0;
    MPI_Allreduce(&localInnerProduct, &globalInnerProduct, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return globalInnerProduct * _params.getEnergyFactor();
}

double SolutionData::getNormL2() const {
    return sqrt(getEnergyL2());
}

double SolutionData::getEnergyL2() const {
    double localEnergy = 0.0;
    const size_t localGridSize = _mpi.getLocalGridSize();

    for (size_t i = 0; i < localGridSize; ++i) {
        localEnergy += norm(_vData[i]);
    }

    double globalEnergy = 0.0;
    MPI_Allreduce(&localEnergy, &globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return globalEnergy * _params.getEnergyFactor();
}

double SolutionData::getEnergyH1() const {
    double localEnergy = 0.0;
    const size_t localGridSize = _mpi.getLocalGridSize();
    vector<double> vH1Weight = _params.getH1Weight();

    for (size_t i = 0; i < localGridSize; ++i) {
        localEnergy += vH1Weight[i] * norm(_vData[i]);
    }
    
    double globalEnergy = 0.0;
    MPI_Allreduce(&localEnergy, &globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return globalEnergy * _params.getEnergyFactor();
}

double SolutionData::getEnergyH2() const {
    double localEnergy = 0.0;
    const size_t localGridSize = _mpi.getLocalGridSize();
    vector<double> vH2Weight = _params.getH2Weight();

    for (size_t i = 0; i < localGridSize; ++i) {
        localEnergy += vH2Weight[i] * norm(_vData[i]);
    }
    
    double globalEnergy = 0.0;
    MPI_Allreduce(&localEnergy, &globalEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return globalEnergy * _params.getEnergyFactor();
}

vector<double> SolutionData::getRadialSpectrum() const {
    vector<double> spectrum(_params.getMaxRadialBin() + 1, 0.0);
    const size_t localGridSize = _mpi.getLocalGridSize();
    vector<size_t> radialBin = _params.getRadialBin();

    for (size_t i = 0; i < localGridSize; ++i) {
        spectrum[radialBin[i]] += norm(_vData[i]);
    }

    MPI_Allreduce(MPI_IN_PLACE, spectrum.data(), static_cast<int>(spectrum.size()), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (double& value : spectrum) {
        value *= 0.5;
    }

    return spectrum;
}

void SolutionData::loadData(SolutionDataType storedDataType, double dCurrentT) {   
    switch (storedDataType) {

        case InitialState: {
            readDistributedFile(_paths.getInitialDataFile());
            break;
        }

        case TerminalState: {
            readDistributedFile(_paths.getTerminalDataFile());
            break;
        }

        case BackwardInitialState: {
            readDistributedFile(_paths.getBackwardSolutionFile());
            break;
        }

        case IntermediateHistory: // fall through  
        
        case RemainderHistory: {
            readDistributedFile(appendTimeStep(_paths.getForwardSolutionFile(), dCurrentT));
            break;
        }
    }
}

void SolutionData::saveData(SolutionDataType storedDataType, double dCurrentT) {   
    switch (storedDataType) {

        case InitialState: {
            writeDistributedFile(_paths.getInitialDataFile());
            break;
        }

        case TerminalState: {
            writeDistributedFile(_paths.getTerminalDataFile());
            break;
        }

        case BackwardInitialState: {
            writeDistributedFile(_paths.getBackwardSolutionFile());
            break;
        }

        case IntermediateHistory: // fall through 

        case RemainderHistory: {
            writeDistributedFile(appendTimeStep(_paths.getForwardSolutionFile(), dCurrentT));
            break;
        }
    }
}

void SolutionData::deleteData() { 
    if (_mpi.isRoot()) {
        filesystem::remove_all(_paths.getDirForwardSolution());
        filesystem::create_directories(_paths.getDirForwardSolution());
    }
    MPI_Barrier(MPI_COMM_WORLD);
}