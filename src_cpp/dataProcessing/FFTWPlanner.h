#ifndef FFTW_H
#define FFTW_H

#include "Parameters.h"

#include <complex>
#include <fftw3.h>

using namespace std;

class FFTWPlanner {
private:
    vector<complex<double>> _dummy;
    fftw_plan _forwardPlan;
    fftw_plan _backwardPlan;

public:
    explicit FFTWPlanner(const Parameters& params);
    ~FFTWPlanner();

    void fft2InPlace(complex<double>* vState);
    void ifft2InPlace(const Parameters& params, complex<double>* vState);
};

#endif