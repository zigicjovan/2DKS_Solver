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
   Nonlinear operator (unchanged, tested)
============================================================ */
static mxArray* compute_nonlinear(
    const mxArray *v_mx,
    const mxArray *D1_mx,
    const mxArray *D2_mx,
    mwSize N)
{
    mxArray *fx_hat_mx = mxCreateDoubleMatrix(N,N,mxCOMPLEX);
    mxArray *fy_hat_mx = mxCreateDoubleMatrix(N,N,mxCOMPLEX);

    const double *v_r  = mxGetPr(v_mx);
    const double *v_i  = mxGetPi(v_mx);
    const double *D1_r = mxGetPr(D1_mx);
    const double *D1_i = mxGetPi(D1_mx);
    const double *D2_r = mxGetPr(D2_mx);
    const double *D2_i = mxGetPi(D2_mx);

    double *fxh_r = mxGetPr(fx_hat_mx);
    double *fxh_i = mxGetPi(fx_hat_mx);
    double *fyh_r = mxGetPr(fy_hat_mx);
    double *fyh_i = mxGetPi(fy_hat_mx);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i = 0; i < C.N2; ++i) {
        double vr=v_r[i], vi=(v_i?v_i[i]:0.0);
        double a1r=D1_r[i], a1i=(D1_i?D1_i[i]:0.0);
        double a2r=D2_r[i], a2i=(D2_i?D2_i[i]:0.0);
        fxh_r[i]=a1r*vr-a1i*vi;
        fxh_i[i]=a1r*vi+a1i*vr;
        fyh_r[i]=a2r*vr-a2i*vi;
        fyh_i[i]=a2r*vi+a2i*vr;
    }

    mxArray *rhs[1], *lhs[1];
    rhs[0]=fx_hat_mx; mexCallMATLAB(1,lhs,1,rhs,"ifft2");
    mxArray *fx_phys=lhs[0];
    rhs[0]=fy_hat_mx; mexCallMATLAB(1,lhs,1,rhs,"ifft2");
    mxArray *fy_phys=lhs[0];

    mxArray *nl_phys=mxCreateDoubleMatrix(N,N,mxREAL);
    double *np=mxGetPr(nl_phys);
    const double *fxp=mxGetPr(fx_phys);
    const double *fyp=mxGetPr(fy_phys);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i=0;i<C.N2;++i)
        np[i]=0.5*(fxp[i]*fxp[i] + fyp[i]*fyp[i]);

    rhs[0]=nl_phys; mexCallMATLAB(1,lhs,1,rhs,"fft2");
    mxArray *nl_hat=lhs[0];

    double *nh_r=mxGetPr(nl_hat);
    double *nh_i=mxGetPi(nl_hat);

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (mwSize i=0;i<C.N2;++i)
        if(!C.mask[i]) { nh_r[i]=0.0; nh_i[i]=0.0; }

    mxDestroyArray(fx_hat_mx);
    mxDestroyArray(fy_hat_mx);
    mxDestroyArray(fx_phys);
    mxDestroyArray(fy_phys);
    mxDestroyArray(nl_phys);

    return nl_hat;
}

/* ============================================================
   MEX ENTRY: FULL RK4 STEP (new interface)
============================================================ */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 10)
        mexErrMsgTxt("Expected 10 inputs.");

    const mxArray *v_mx   = prhs[0];
    const mxArray *D1_mx  = prhs[1];
    const mxArray *D2_mx  = prhs[2];
    mwSize N              = (mwSize)mxGetScalar(prhs[3]);
    const mxArray *Lin_mx = prhs[4];
    double dt             = mxGetScalar(prhs[5]);
    const double *alpha_I = mxGetPr(prhs[6]);
    const double *beta_I  = mxGetPr(prhs[7]);
    const double *alpha_E = mxGetPr(prhs[8]);
    const double *beta_E  = mxGetPr(prhs[9]);

    ensure_cache(N);

    const double *Lin_r = mxGetPr(Lin_mx);
    const double *Lin_i = mxGetPi(Lin_mx);  // may be NULL if Lin is real

    mxArray *v_step = mxDuplicateArray(v_mx);

    /* initialize Nonlin_v0 internally */
    mxArray *Nonlin_v0 = compute_nonlinear(v_step,D1_mx,D2_mx,N);

    for (int k=0;k<4;++k) {

        mxArray *Nonlin_v1 = compute_nonlinear(v_step,D1_mx,D2_mx,N);

        double *v_r = mxGetPr(v_step);
        double *v_i = mxGetPi(v_step);
        const double *n1_r = mxGetPr(Nonlin_v1);
        const double *n1_i = mxGetPi(Nonlin_v1);
        const double *n0_r = mxGetPr(Nonlin_v0);
        const double *n0_i = mxGetPi(Nonlin_v0);

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (mwSize i=0;i<C.N2;++i) {
            double Lr = Lin_r[i];
            double Li = (Lin_i ? Lin_i[i] : 0.0);
            
            double vr = v_r[i];
            double vi = (v_i ? v_i[i] : 0.0);
            
            /* A = 1 - dt*beta_I[k]*L */
            double Ar_r = 1.0 - dt*beta_I[k]*Lr;
            double Ar_i =       - dt*beta_I[k]*Li;
            
            /* B = 1 + dt*alpha_I[k]*L */
            double Br_r = 1.0 + dt*alpha_I[k]*Lr;
            double Br_i =       dt*alpha_I[k]*Li;
            
            /* A*v (complex multiply) */
            double Av_r = Ar_r*vr - Ar_i*vi;
            double Av_i = Ar_r*vi + Ar_i*vr;
            
            /* rhs = A*v - dt*alpha_E*n1 - dt*beta_E*n0 */
            double rhsr = Av_r - dt*alpha_E[k]*n1_r[i] - dt*beta_E[k]*n0_r[i];
            double rhsi = Av_i - dt*alpha_E[k]*n1_i[i] - dt*beta_E[k]*n0_i[i];
            
            /* v_new = rhs / B (complex divide) */
            double denom = Br_r*Br_r + Br_i*Br_i;
            double vnr = (rhsr*Br_r + rhsi*Br_i) / denom;
            double vni = (rhsi*Br_r - rhsr*Br_i) / denom;
            
            v_r[i] = vnr;
            if (v_i) v_i[i] = vni;
        }

        mxDestroyArray(Nonlin_v0);
        Nonlin_v0 = Nonlin_v1;
    }

    plhs[0]=v_step;
}
