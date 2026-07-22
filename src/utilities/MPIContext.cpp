#include "MPIContext.h"

#include <fftw3-mpi.h>
#include <iostream>

using namespace std;

MPIContext::MPIContext(std::size_t gridSize1, std::size_t gridSize2) {
    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_size);

    _gridSize1 = gridSize1;
    _localAllocationSize = fftw_mpi_local_size_2d(static_cast<std::ptrdiff_t>(gridSize2), static_cast<std::ptrdiff_t>(gridSize1),
                                                  MPI_COMM_WORLD, &_localGridSize2, &_localGridStart2);
}

int MPIContext::getRank() const {
    return _rank;
}

int MPIContext::getSize() const {
    return _size;
}

std::size_t MPIContext::getLocalGridSize2() const {
    return static_cast<std::size_t>(_localGridSize2);
}

std::size_t MPIContext::getLocalGridStart2() const {
    return static_cast<std::size_t>(_localGridStart2);
}

std::size_t MPIContext::getLocalGridSize() const {
    return static_cast<std::size_t>(_localGridSize2) * _gridSize1;
}

std::size_t MPIContext::getLocalAllocationSize() const {
    return static_cast<std::size_t>(_localAllocationSize);
}

void MPIContext::printDecomposition() const {
    for (int rank = 0; rank < _size; ++rank) {
        if (_rank == rank) {
            const size_t rowStart = getLocalGridStart2();
            const size_t rowEnd = rowStart + getLocalGridSize2();
            cout << "Rank " << _rank << " of " << _size << ": local rows = " << getLocalGridSize2() << ", global row range = [" << rowStart
                      << ", " << rowEnd << "), local grid elements = " << getLocalGridSize() << '\n' << flush;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

bool MPIContext::isRoot() const {
    return _rank == 0;
}