#ifndef FFTW_H
#define FFTW_H

#include "Parameters.h"

#include <complex>
#include <fftw3.h>

using Complex = std::complex<double>;

class FFTWPlanner {
private:
    std::vector<Complex> dummy;
public:
    fftw_plan forwardPlan;
    fftw_plan backwardPlan;
    FFTWPlanner(const Parameters& params) {
        dummy.resize(params.iTotalGridSize);

        forwardPlan = fftw_plan_dft_2d(
            params.iGridSize2, params.iGridSize1,
            reinterpret_cast<fftw_complex*>(dummy.data()),
            reinterpret_cast<fftw_complex*>(dummy.data()),
            FFTW_FORWARD, FFTW_MEASURE);

        backwardPlan = fftw_plan_dft_2d(
            params.iGridSize2, params.iGridSize1,
            reinterpret_cast<fftw_complex*>(dummy.data()),
            reinterpret_cast<fftw_complex*>(dummy.data()),
            FFTW_BACKWARD, FFTW_MEASURE);
    }
    ~FFTWPlanner() {
        fftw_destroy_plan(forwardPlan);
        fftw_destroy_plan(backwardPlan);
    }

    void fft2InPlace(Complex* vState);
    void ifft2InPlace(const Parameters& params, Complex* vState);
};

#endif