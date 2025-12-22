#include "mex.h"
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

/* ============================================================
   Cache: dealias mask only
============================================================ */
struct Cache {
    mwSize N  = 0;
    mwSize N2 = 0;
    std::vector<char> mask;
};

static Cache C;

/* ============================================================
   Build 2/3 dealias mask
============================================================ */
static void build_mask(mwSize N)
{
    C.mask.assign(N*N, 0);
    double kcut = (2.0/3.0) * (N/2.0);

    std::vector<double> k(N);
    for (mwSize i = 0; i < N; ++i)
        k[i] = (i <= N/2) ? i : (double)i - (double)N;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize j = 0; j < N; ++j)
        for (mwSize i = 0; i < N; ++i)
            if (std::abs(k[i]) <= kcut && std::abs(k[j]) <= kcut)
                C.mask[i + j*N] = 1;
}

static void ensure_cache(mwSize N)
{
    if (C.N == N) return;
    C.N = N;
    C.N2 = N*N;
    build_mask(N);
}

/* ============================================================
   MEX ENTRY
============================================================ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgTxt("Usage: ks_nonlinear_mex(v,D1,D2,N)");

    mwSize N = (mwSize)mxGetScalar(prhs[3]);
    ensure_cache(N);

    const mxArray *v_mx  = prhs[0];
    const mxArray *D1_mx = prhs[1];
    const mxArray *D2_mx = prhs[2];

    if (mxGetNumberOfElements(v_mx) != N*N)
        mexErrMsgTxt("v must be N^2 elements.");

#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *v  = mxGetComplexDoubles(v_mx);
    const mxComplexDouble *D1 = mxGetComplexDoubles(D1_mx);
    const mxComplexDouble *D2 = mxGetComplexDoubles(D2_mx);
#else
    const double *v_r  = mxGetPr(v_mx);
    const double *v_i  = mxGetPi(v_mx);
    const double *D1_r = mxGetPr(D1_mx);
    const double *D1_i = mxGetPi(D1_mx);
    const double *D2_r = mxGetPr(D2_mx);
    const double *D2_i = mxGetPi(D2_mx);
#endif

    /* ------------------------------------------------------------
       1) Spectral derivatives
    ------------------------------------------------------------ */
    mxArray *fx_hat_mx = mxCreateDoubleMatrix(N,N,mxCOMPLEX);
    mxArray *fy_hat_mx = mxCreateDoubleMatrix(N,N,mxCOMPLEX);

#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *fxh = mxGetComplexDoubles(fx_hat_mx);
    mxComplexDouble *fyh = mxGetComplexDoubles(fy_hat_mx);
#else
    double *fxh_r = mxGetPr(fx_hat_mx);
    double *fxh_i = mxGetPi(fx_hat_mx);
    double *fyh_r = mxGetPr(fy_hat_mx);
    double *fyh_i = mxGetPi(fy_hat_mx);
#endif

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i = 0; i < C.N2; ++i) {

#if MX_HAS_INTERLEAVED_COMPLEX
        double vr = v[i].real, vi = v[i].imag;
        double a1r = D1[i].real, a1i = D1[i].imag;
        double a2r = D2[i].real, a2i = D2[i].imag;

        fxh[i].real = a1r*vr - a1i*vi;
        fxh[i].imag = a1r*vi + a1i*vr;
        fyh[i].real = a2r*vr - a2i*vi;
        fyh[i].imag = a2r*vi + a2i*vr;
#else
        double vr  = v_r[i];
        double vi  = v_i ? v_i[i] : 0.0;
        double a1r = D1_r[i];
        double a1i = D1_i ? D1_i[i] : 0.0;
        double a2r = D2_r[i];
        double a2i = D2_i ? D2_i[i] : 0.0;

        fxh_r[i] = a1r*vr - a1i*vi;
        fxh_i[i] = a1r*vi + a1i*vr;
        fyh_r[i] = a2r*vr - a2i*vi;
        fyh_i[i] = a2r*vi + a2i*vr;
#endif
    }

    /* ------------------------------------------------------------
       2) Inverse FFTs (MATLAB)
    ------------------------------------------------------------ */
    mxArray *rhs[1], *lhs[1];

    rhs[0] = fx_hat_mx;
    mexCallMATLAB(1, lhs, 1, rhs, "ifft2");
    mxArray *fx_phys_mx = lhs[0];

    rhs[0] = fy_hat_mx;
    mexCallMATLAB(1, lhs, 1, rhs, "ifft2");
    mxArray *fy_phys_mx = lhs[0];

#if MX_HAS_INTERLEAVED_COMPLEX
    const mxComplexDouble *fxp = mxGetComplexDoubles(fx_phys_mx);
    const mxComplexDouble *fyp = mxGetComplexDoubles(fy_phys_mx);
#else
    const double *fxp_r = mxGetPr(fx_phys_mx);
    const double *fyp_r = mxGetPr(fy_phys_mx);
#endif

    /* ------------------------------------------------------------
       3) Nonlinear term
    ------------------------------------------------------------ */
    mxArray *nonlin_phys_mx = mxCreateDoubleMatrix(N,N,mxREAL);
    double *np = mxGetPr(nonlin_phys_mx);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i = 0; i < C.N2; ++i) {
#if MX_HAS_INTERLEAVED_COMPLEX
        double fx = fxp[i].real;
        double fy = fyp[i].real;
#else
        double fx = fxp_r[i];
        double fy = fyp_r[i];
#endif
        np[i] = 0.5*(fx*fx + fy*fy);
    }

    /* ------------------------------------------------------------
       4) Forward FFT
    ------------------------------------------------------------ */
    rhs[0] = nonlin_phys_mx;
    mexCallMATLAB(1, lhs, 1, rhs, "fft2");
    mxArray *nonlin_hat_mx = lhs[0];

#if MX_HAS_INTERLEAVED_COMPLEX
    mxComplexDouble *nh = mxGetComplexDoubles(nonlin_hat_mx);
#else
    double *nh_r = mxGetPr(nonlin_hat_mx);
    double *nh_i = mxGetPi(nonlin_hat_mx);
#endif

    /* ------------------------------------------------------------
       5) Dealias
    ------------------------------------------------------------ */
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i = 0; i < C.N2; ++i)
        if (!C.mask[i]) {
#if MX_HAS_INTERLEAVED_COMPLEX
            nh[i].real = 0.0;
            nh[i].imag = 0.0;
#else
            nh_r[i] = 0.0;
            nh_i[i] = 0.0;
#endif
        }

    /* ------------------------------------------------------------
       Output
    ------------------------------------------------------------ */
    plhs[0] = nonlin_hat_mx;

    mxDestroyArray(fx_hat_mx);
    mxDestroyArray(fy_hat_mx);
    mxDestroyArray(fx_phys_mx);
    mxDestroyArray(fy_phys_mx);
    mxDestroyArray(nonlin_phys_mx);
}
