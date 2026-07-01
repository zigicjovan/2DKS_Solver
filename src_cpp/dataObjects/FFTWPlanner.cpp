#include "Parameters.h"
#include "FFTWPlanner.h"
#include <complex>
#include <fftw3.h>

using Complex = std::complex<double>;

void FFTWPlanner::fft2InPlace(Complex* vState) {
    fftw_execute_dft(forwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));
}

void FFTWPlanner::ifft2InPlace(const Parameters& params, Complex* vState) {
    fftw_execute_dft(backwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));

    const double scale = 1.0 / static_cast<double>(params.iTotalGridSize);
    for (std::size_t p = 0; p < params.iTotalGridSize; ++p)
        vState[p] *= scale;
}