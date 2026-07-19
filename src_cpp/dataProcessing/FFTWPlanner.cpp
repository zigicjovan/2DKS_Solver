#include "Parameters.h"
#include "FFTWPlanner.h"
#include <complex>
#include <fftw3.h>

using namespace std;

FFTWPlanner::FFTWPlanner(const Parameters& params) {
    _dummy.resize(params.getTotalGridSize());

    _forwardPlan = fftw_plan_dft_2d( params.getGridSize2(), params.getGridSize1(),
        reinterpret_cast<fftw_complex*>(_dummy.data()), reinterpret_cast<fftw_complex*>(_dummy.data()), FFTW_FORWARD, FFTW_MEASURE);

    _backwardPlan = fftw_plan_dft_2d( params.getGridSize2(), params.getGridSize1(),
        reinterpret_cast<fftw_complex*>(_dummy.data()), reinterpret_cast<fftw_complex*>(_dummy.data()), FFTW_BACKWARD, FFTW_MEASURE);
}

FFTWPlanner::~FFTWPlanner() {
    fftw_destroy_plan(_forwardPlan);
    fftw_destroy_plan(_backwardPlan);
}

void FFTWPlanner::fft2InPlace(complex<double>* vState) {
    fftw_execute_dft(_forwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));
}

void FFTWPlanner::ifft2InPlace(const Parameters& params, complex<double>* vState) {
    fftw_execute_dft(_backwardPlan,reinterpret_cast<fftw_complex*>(vState),reinterpret_cast<fftw_complex*>(vState));

    const double scale = 1.0 / static_cast<double>(params.getTotalGridSize());
    size_t gridSize = params.getTotalGridSize();
    
    for (size_t i = 0; i < gridSize; ++i) {
        vState[i] *= scale;
    }
}