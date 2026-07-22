#ifndef FFTW_H
#define FFTW_H

#include "Parameters.h"
#include "MPIContext.h"

#include <fftw3-mpi.h>
#include <complex>

using namespace std;

class FFTWPlanner {
private:
    const MPIContext& _mpi;

    vector<complex<double>> _dummy;
    fftw_plan _forwardPlan;
    fftw_plan _backwardPlan;

public:
    FFTWPlanner(const Parameters& params, const MPIContext& mpi);
    ~FFTWPlanner();

    void fft2InPlace(complex<double>* vState);
    void ifft2InPlace(const Parameters& params, complex<double>* vState);
};

#endif