#include "Parameters.h"
#include "FFTWPlanner.h"
#include <complex>
#include <fftw3.h>

using namespace std;

void FFTWPlanner::fft2InPlace(complex<double>* vState) {
    fftw_execute_dft(forwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));
}

void FFTWPlanner::ifft2InPlace(const Parameters& params, complex<double>* vState) {
    fftw_execute_dft(backwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));

    const double scale = 1.0 / static_cast<double>(params.iTotalGridSize);
    for (size_t i = 0; i < params.iTotalGridSize; ++i)
        vState[i] *= scale;
}