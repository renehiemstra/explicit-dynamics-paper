/*==========================================================
 * quadratureLoop.c - MATLAB External Interface
*
* Performs the quadrature loop to fill objects for use in 
* sum factorization, accompanying the provided matlab code.
*
* The calling syntax is:
*
*		void quadratureLoop(m1, m2, m3
                            sumfact,
                            weights,
                            jacobian)
*
* This is a MEX-file for MATLAB.
* Copyright 2007-2012 The MathWorks, Inc.
*
*========================================================*/

#include "mex.h"
#include "math.h"

double det_mat_2d(double a, double b, double c, double d){
    return a * d - b * c;
}

double det_mat_3d(double A[3][3]){
    // compute transpose of the cofactor matrix
    double d = 0.0;
    double s = 1;
    for (mwSize i=0; i<3; i++) {
        d = d + s * A[i][0] * det_mat_2d(A[(i+1)%3][1], A[(i+1)%3][2], A[(i+2)%3][1], A[(i+2)%3][2]);
        s = -s;
    }
    return d;
}

double norm(double a1, double a2, double a3){
    return sqrt(a1*a1 + a2*a2 + a3*a3);
}

/* The computational routine */
void quadrature_loop(mwSize m1, mwSize m2, mwSize m3,
            double *sumfact[4],
            double *weights[3],
            double *jacobian[3][3])
{
    mwSize k;
    mwSize l;
    mwSize m;
    mwSize lindex;

    double w, c;
    double A[3][3];

    // loop over points in w-dir
    for (m=0; m<m3; m++) {
        // loop over points in v-dir
        for (l=0; l<m2; l++) {
            // loop over points in u-dir
            for (k=0; k<m1; k++){
                // compute linear index
                lindex = m1 * m2 * m + m1 * l + k;
    
                // load precomputed jacobian matrix terms
                A[0][0] = jacobian[0][0][lindex]; A[0][1] = jacobian[0][1][lindex]; A[0][2] = jacobian[0][2][lindex];
                A[1][0] = jacobian[1][0][lindex]; A[1][1] = jacobian[1][1][lindex]; A[1][2] = jacobian[1][2][lindex];
                A[2][0] = jacobian[2][0][lindex]; A[2][1] = jacobian[2][1][lindex]; A[2][2] = jacobian[2][2][lindex];
                c = det_mat_3d(A);

                // load quadrature weight
                w = weights[0][k] * weights[1][l] * weights[2][m] * c;
    
                // add results to sumfactorization objects.
                // the terms in g are weighted by the quaderature weight.
                sumfact[0][lindex] = w;
                
            } // loop over 'k'
        } // loop over 'l'
    } // loop over 'm'
} // end quadratureLoop


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // sizes
    mwSize m1;
    mwSize m2;
    mwSize m3;
    
    // pointer objects
    double *sumfact[4];
    double *weights[3];
    double *jacobian[3][3];

    /* get the value of the scalar input  */
    m1 = mxGetScalar(prhs[0]);
    m2 = mxGetScalar(prhs[1]);
    m3 = mxGetScalar(prhs[2]);
    
    /* create a pointer to the real data in the input matrix  */
    sumfact[0] = mxGetDoubles(mxGetCell(prhs[3], 0)); 
    sumfact[1] = mxGetDoubles(mxGetCell(prhs[3], 1)); 
    sumfact[2] = mxGetDoubles(mxGetCell(prhs[3], 2));
    sumfact[3] = mxGetDoubles(mxGetCell(prhs[3], 3));

    weights[0] = mxGetDoubles(mxGetCell(prhs[4], 0)); 
    weights[1] = mxGetDoubles(mxGetCell(prhs[4], 1));
    weights[2] = mxGetDoubles(mxGetCell(prhs[4], 2));

    jacobian[0][0] = mxGetDoubles(mxGetCell(prhs[5], 0)); 
    jacobian[1][0] = mxGetDoubles(mxGetCell(prhs[5], 1));
    jacobian[2][0] = mxGetDoubles(mxGetCell(prhs[5], 2));
    jacobian[0][1] = mxGetDoubles(mxGetCell(prhs[5], 3)); 
    jacobian[1][1] = mxGetDoubles(mxGetCell(prhs[5], 4));
    jacobian[2][1] = mxGetDoubles(mxGetCell(prhs[5], 5));
    jacobian[0][2] = mxGetDoubles(mxGetCell(prhs[5], 6)); 
    jacobian[1][2] = mxGetDoubles(mxGetCell(prhs[5], 7));
    jacobian[2][2] = mxGetDoubles(mxGetCell(prhs[5], 8));

    /* call the computational routine */
    quadrature_loop(m1, m2, m3,
        sumfact, weights,
                jacobian);
}
