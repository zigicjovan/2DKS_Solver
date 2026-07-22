#include "Parameters.h"
#include "FFTWPlanner.h"

#include <complex>
#include <fftw3-mpi.h>

using namespace std;

FFTWPlanner::FFTWPlanner(const Parameters& params, const MPIContext& mpi) : _mpi(mpi) {
    _dummy.resize(_mpi.getLocalAllocationSize());

    _forwardPlan = fftw_mpi_plan_dft_2d(static_cast<ptrdiff_t>(params.getGridSize2()), static_cast<ptrdiff_t>(params.getGridSize1()),
                                        reinterpret_cast<fftw_complex*>(_dummy.data()), reinterpret_cast<fftw_complex*>(_dummy.data()),
                                        MPI_COMM_WORLD, FFTW_FORWARD, FFTW_MEASURE);

    _backwardPlan = fftw_mpi_plan_dft_2d(static_cast<ptrdiff_t>(params.getGridSize2()), static_cast<ptrdiff_t>(params.getGridSize1()),
                                         reinterpret_cast<fftw_complex*>(_dummy.data()), reinterpret_cast<fftw_complex*>(_dummy.data()),
                                         MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_MEASURE);
}

FFTWPlanner::~FFTWPlanner() {
    fftw_destroy_plan(_forwardPlan);
    fftw_destroy_plan(_backwardPlan);
}

void FFTWPlanner::fft2InPlace(complex<double>* vState) {
    fftw_mpi_execute_dft(_forwardPlan, reinterpret_cast<fftw_complex*>(vState), reinterpret_cast<fftw_complex*>(vState));
}

void FFTWPlanner::ifft2InPlace(const Parameters& params, complex<double>* vState) {
    fftw_mpi_execute_dft(_backwardPlan, reinterpret_cast<fftw_complex*>(vState), reinterpret_cast<fftw_complex*>(vState));

    const double scale = 1.0 / static_cast<double>(params.getTotalGridSize());
    const size_t localGridSize = _mpi.getLocalGridSize();
    
    for (size_t i = 0; i < localGridSize; ++i) {
        vState[i] *= scale;
    }
}