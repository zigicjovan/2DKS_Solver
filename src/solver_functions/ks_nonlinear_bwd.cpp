#include "mex.h"
#include <vector>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

/* =========================================================
   Cache
========================================================= */
struct Cache {
    mwSize N = 0, N2 = 0;
    std::vector<char> mask;
};
static Cache C;

static void build_mask(mwSize N) {
    C.mask.resize(N*N);
    double kcut = (2.0/3.0)*(N/2.0);

    std::vector<double> k(N);
    for (mwSize i=0;i<N;i++)
        k[i] = (i<=N/2)? (double)i : (double)i - (double)N;

    #pragma omp parallel for
    for (mwSize j=0;j<N;j++) {
        for (mwSize i=0;i<N;i++) {
            C.mask[i + j*N] =
                (std::abs(k[i])<=kcut && std::abs(k[j])<=kcut);
        }
    }
}

static void ensure_cache(mwSize N) {
    if (C.N == N) return;
    C = Cache();
    C.N = N;
    C.N2 = N*N;
    build_mask(N);
}

/* =========================================================
   Helpers: hard dimension checks
========================================================= */
static void require_numel(const mxArray* A, mwSize N2, const char* name) {
    if (mxGetNumberOfElements(A) != N2) {
        mexErrMsgIdAndTxt("ks:dim",
            "%s must have numel == N^2. Got %llu, expected %llu.",
            name,
            (unsigned long long)mxGetNumberOfElements(A),
            (unsigned long long)N2);
    }
}

static void require_real_double_vector4(const mxArray* A, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A) || mxGetNumberOfElements(A) != 4) {
        mexErrMsgIdAndTxt("ks:dim",
            "%s must be a real double vector with 4 elements.", name);
    }
}

static void require_real_double_lenN2(const mxArray* A, mwSize N2, const char* name) {
    if (!mxIsDouble(A) || mxIsComplex(A) || mxGetNumberOfElements(A) != N2) {
        mexErrMsgIdAndTxt("ks:dim",
            "%s must be a real double array with numel == N^2.", name);
    }
}

/* =========================================================
   MEX ENTRY
========================================================= */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    if (nrhs != 12) {
        mexErrMsgTxt("Usage: v1 = ks_nonlinear_bwd(z, fwd, D1, D2, Lap, N, Lin, dt, aI, bI, aE, bE)");
    }

    const mxArray *z_hat = prhs[0];     // complex, N^2
    const mxArray *fwd   = prhs[1];     // complex, N^2
    const mxArray *D1    = prhs[2];     // complex, N^2
    const mxArray *D2    = prhs[3];     // complex, N^2
    const mxArray *Lap   = prhs[4];     // complex, N^2
    mwSize N             = (mwSize)mxGetScalar(prhs[5]);
    const mxArray *Lin   = prhs[6];     // real,   N^2
    double dt            = mxGetScalar(prhs[7]);
    const mxArray *aI_mx = prhs[8];
    const mxArray *bI_mx = prhs[9];
    const mxArray *aE_mx = prhs[10];
    const mxArray *bE_mx = prhs[11];

    ensure_cache(N);
    mwSize N2 = C.N2;

    /* ---- hard dimension checks ---- */
    require_numel(z_hat, N2, "z");
    require_numel(fwd,   N2, "fwd");
    require_numel(D1,    N2, "D1");
    require_numel(D2,    N2, "D2");
    require_numel(Lap,   N2, "Lap");
    require_real_double_lenN2(Lin, N2, "Lin");
    require_real_double_vector4(aI_mx, "alpha_I");
    require_real_double_vector4(bI_mx, "beta_I");
    require_real_double_vector4(aE_mx, "alpha_E");
    require_real_double_vector4(bE_mx, "beta_E");

    const double *aI = mxGetPr(aI_mx);
    const double *bI = mxGetPr(bI_mx);
    const double *aE = mxGetPr(aE_mx);
    const double *bE = mxGetPr(bE_mx);
    const double *Linr = mxGetPr(Lin);

    /* ---- allocate stage variables (vector form: N^2 x 1) ---- */
    mxArray *zcur = mxDuplicateArray(z_hat);                 // complex
    mxArray *NL0  = mxCreateDoubleMatrix(N2,1,mxCOMPLEX);    // complex zeros

    /* ---- precompute forward derivatives in Fourier: D1.*fwd, D2.*fwd, Lap.*fwd ---- */
    mxArray *fwd_fx_hat  = mxDuplicateArray(fwd);
    mxArray *fwd_fy_hat  = mxDuplicateArray(fwd);
    mxArray *fwd_lap_hat = mxDuplicateArray(fwd);

    double *fxr = mxGetPr(fwd_fx_hat),  *fxi = mxGetPi(fwd_fx_hat);
    double *fyr = mxGetPr(fwd_fy_hat),  *fyi = mxGetPi(fwd_fy_hat);
    double *flr = mxGetPr(fwd_lap_hat), *fli = mxGetPi(fwd_lap_hat);

    const double *D1r = mxGetPr(D1),  *D1i = mxGetPi(D1);
    const double *D2r = mxGetPr(D2),  *D2i = mxGetPi(D2);
    const double *Lr  = mxGetPr(Lap), *Li  = mxGetPi(Lap);

    #pragma omp parallel for
    for (mwSize i=0;i<N2;i++) {
        /* fx_hat = D1 .* fwd */
        double vr = fxr[i], vi = (fxi?fxi[i]:0.0);
        double a1r = D1r[i], a1i = (D1i?D1i[i]:0.0);
        fxr[i] = a1r*vr - a1i*vi;
        if (fxi) fxi[i] = a1r*vi + a1i*vr;

        /* fy_hat = D2 .* fwd */
        vr = fyr[i]; vi = (fyi?fyi[i]:0.0);
        double a2r = D2r[i], a2i = (D2i?D2i[i]:0.0);
        fyr[i] = a2r*vr - a2i*vi;
        if (fyi) fyi[i] = a2r*vi + a2i*vr;

        /* lap_hat = Lap .* fwd */
        vr = flr[i]; vi = (fli?fli[i]:0.0);
        double alr = Lr[i], ali = (Li?Li[i]:0.0);
        flr[i] = alr*vr - ali*vi;
        if (fli) fli[i] = alr*vi + ali*vr;
    }

    /* =========================================================
       RK4 loop
    ========================================================= */
    for (int k=0;k<4;k++) {

        /* z_x_hat = D1 .* zcur,  z_y_hat = D2 .* zcur */
        mxArray *zx_hat = mxDuplicateArray(zcur);
        mxArray *zy_hat = mxDuplicateArray(zcur);

        double *zxr = mxGetPr(zx_hat), *zxi = mxGetPi(zx_hat);
        double *zyr = mxGetPr(zy_hat), *zyi = mxGetPi(zy_hat);

        double *zhr = mxGetPr(zcur);
        double *zhi = mxGetPi(zcur);

        #pragma omp parallel for
        for (mwSize i=0;i<N2;i++) {
            double vr = zhr[i];
            double vi = (zhi?zhi[i]:0.0);

            double a1r = D1r[i], a1i = (D1i?D1i[i]:0.0);
            zxr[i] = a1r*vr - a1i*vi;
            if (zxi) zxi[i] = a1r*vi + a1i*vr;

            double a2r = D2r[i], a2i = (D2i?D2i[i]:0.0);
            zyr[i] = a2r*vr - a2i*vi;
            if (zyi) zyi[i] = a2r*vi + a2i*vr;
        }

        /* ---- physical-space fields: real(ifft2(.)) ---- */
        mxArray *zx_phys = nullptr, *zy_phys = nullptr;
        mxArray *fx_phys = nullptr, *fy_phys = nullptr;
        mxArray *lp_phys = nullptr, *z_phys  = nullptr;

        mexCallMATLAB(1,&zx_phys,1,&zx_hat,      "ifft2");
        mexCallMATLAB(1,&zy_phys,1,&zy_hat,      "ifft2");
        mexCallMATLAB(1,&fx_phys,1,&fwd_fx_hat,  "ifft2");
        mexCallMATLAB(1,&fy_phys,1,&fwd_fy_hat,  "ifft2");
        mexCallMATLAB(1,&lp_phys,1,&fwd_lap_hat, "ifft2");
        mexCallMATLAB(1,&z_phys, 1,&zcur,        "ifft2");

        /* Use real parts only (consistent with your MATLAB multiply2D 'fourier2real') */
        const double *zxP = mxGetPr(zx_phys);
        const double *zyP = mxGetPr(zy_phys);
        const double *fxP = mxGetPr(fx_phys);
        const double *fyP = mxGetPr(fy_phys);
        const double *lpP = mxGetPr(lp_phys);
        const double *zP  = mxGetPr(z_phys);

        /* ---- nonlinear term in physical space (real N x N) ---- */
        mxArray *NL_phys = mxCreateDoubleMatrix(N, N, mxREAL);
        double *nl = mxGetPr(NL_phys);

        #pragma omp parallel for
        for (mwSize i=0;i<N2;i++) {
            nl[i] = -( zxP[i]*fxP[i] + zyP[i]*fyP[i] ) - ( lpP[i]*zP[i] );
        }

        /* ---- FFT + dealias ---- */
        mxArray *NL_hat = nullptr;
        mexCallMATLAB(1,&NL_hat,1,&NL_phys,"fft2");

        double *nr = mxGetPr(NL_hat);
        double *ni = mxGetPi(NL_hat);

        #pragma omp parallel for
        for (mwSize i=0;i<N2;i++) {
            if (!C.mask[i]) {
                nr[i] = 0.0;
                if (ni) ni[i] = 0.0;
            }
        }

        /* ---- IMEX update in Fourier: divide by (1 + dt*alpha_I(k)*Lin) ---- */
        double *n0r = mxGetPr(NL0);
        double *n0i = mxGetPi(NL0);
        double *n1r = mxGetPr(NL_hat);
        double *n1i = mxGetPi(NL_hat);

        #pragma omp parallel for
        for (mwSize i=0;i<N2;i++) {
            double L = Linr[i];
            double Ar = 1.0 - dt*bI[k]*L;
            double Br = 1.0 + dt*aI[k]*L;

            double vr = zhr[i];
            double vi = (zhi?zhi[i]:0.0);

            double rhsr = Ar*vr - dt*aE[k]*n1r[i] - dt*bE[k]*n0r[i];
            double rhsi = Ar*vi - dt*aE[k]*(n1i?n1i[i]:0.0) - dt*bE[k]*(n0i?n0i[i]:0.0);

            /* divide complex by real Br */
            zhr[i] = rhsr / Br;
            if (zhi) zhi[i] = rhsi / Br;
        }

        /* stage carry */
        mxDestroyArray(NL0);
        NL0 = NL_hat;

        /* cleanup temporaries */
        mxDestroyArray(zx_hat);
        mxDestroyArray(zy_hat);
        mxDestroyArray(zx_phys);
        mxDestroyArray(zy_phys);
        mxDestroyArray(fx_phys);
        mxDestroyArray(fy_phys);
        mxDestroyArray(lp_phys);
        mxDestroyArray(z_phys);
        mxDestroyArray(NL_phys);
    }

    /* output */
    plhs[0] = zcur;

    /* cleanup persistent locals */
    mxDestroyArray(NL0);
    mxDestroyArray(fwd_fx_hat);
    mxDestroyArray(fwd_fy_hat);
    mxDestroyArray(fwd_lap_hat);
}
