#ifndef MPICONTEXT_H
#define MPICONTEXT_H

#include <cstddef>
#include <mpi.h>

class MPIContext {
private:
    int _rank;
    int _size;

    std::size_t _gridSize1;
    std::ptrdiff_t _localGridSize2;
    std::ptrdiff_t _localGridStart2;
    std::ptrdiff_t _localAllocationSize;

public:
    MPIContext(std::size_t gridSize1, std::size_t gridSize2);

    int getRank() const;
    int getSize() const;

    std::size_t getLocalGridSize2() const;
    std::size_t getLocalGridStart2() const;
    std::size_t getLocalGridSize() const;
    std::size_t getLocalAllocationSize() const;
    void printDecomposition() const;

    bool isRoot() const;
};

#endif