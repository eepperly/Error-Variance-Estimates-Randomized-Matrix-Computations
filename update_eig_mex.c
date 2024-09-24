// Compute eigendecomposition of diag(d) + c*v*v', for a unit vector v

#include "mex.h"
#include "lapack.h"
#include "math.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#if MX_HAS_INTERLEAVED_COMPLEX
    mxDouble *d, *v, *c, *Q, *l;
#else
    double *d, *v, *c, *Q, *l;
#endif
    ptrdiff_t n;      /* matrix dimensions */

#if MX_HAS_INTERLEAVED_COMPLEX
    d = mxGetDoubles(prhs[0]);
    v = mxGetDoubles(prhs[1]);
    c = mxGetDoubles(prhs[2]);
#else
    d = mxGetPr(prhs[0]);
    v = mxGetPr(prhs[1]);
    c = mxGetPr(prhs[2]);
#endif
    n = mxGetM(prhs[0]);  

    for (int i = 0; i < 3; ++i) {
        if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) {
          mexErrMsgIdAndTxt( "MATLAB:update_eig:fieldNotRealMatrix",
                  "Input arguments must be a real matrix.");
        }
    }

    if (mxGetN(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:update_eig:dims",
                "First input must be column vector");
    }

    if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:update_eig:dims",
                "Second input must be column vector, same size as first");
    }


    if (mxGetM(prhs[2]) != 1 || mxGetN(prhs[1]) != 1) {
        mexErrMsgIdAndTxt("MATLAB:update_eig:dims",
                "Third input must be scalar");
    }

    /* create output matrix C */
    plhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);;
#if MX_HAS_INTERLEAVED_COMPLEX
    Q = mxGetDoubles(plhs[0]);
    l = mxGetDoubles(plhs[1]);
#else
    Q = mxGetPr(plhs[0]);
    l = mxGetPr(plhs[1]);
#endif

    double *work = malloc(sizeof(double) * n*n);

    ptrdiff_t info;
    ptrdiff_t one = 1;

    /* Pass arguments to Fortran by reference */
    dlaed9(&n, &one, &n, &n, l, work, &n, c, d, v, Q, &n, &info);

    if (info != 0) {
        if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != 1) {
            mexWarnMsgIdAndTxt("MATLAB:update_eig:dslaed9_error",
                    "dlaed9 exited with error %d", info);
        }
    }

    bool found_nans = false;
    for (int i = 0; i < n; ++i){
        if (isnan(l[i])) {
            mexWarnMsgIdAndTxt("MATLAB:update_eig:nans",
                "d vector has NaNs, problem not properly deflated");
            break;
        }
        for (int j = 0; j < n; ++j) {
            if (isnan(Q[i + n*j])) {
                mexWarnMsgIdAndTxt("MATLAB:update_eig:nans",
                    "Q matrix has NaNs, problem not properly deflated");
                found_nans = true;
                break;
            }
        }
        if (found_nans) break;
    }

    free(work);
}
