/*==========================================================
 * quadratureLoop.c - MATLAB External Interface
*
* Performs the quadrature loop to fill objects for use in 
* sum factorization, accompanying the provided matlab code.
*
* The calling syntax is:
*
*		void quadratureLoop(m1, m2,
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

/* The computational routine */
void quadrature_loop(mwSize m1, mwSize m2,
            double *sumfact[3],
            double *weights[2],
            double *jacobian[2][2])
{
    mwSize k;
    mwSize l;
    mwSize lindex;

    double w, c;
    
    // loop over points in v-dir
    for (l=0; l<m2; l++) {
        // loop over points in u-dir
        for (k=0; k<m1; k++){
            // compute linear index
            lindex = m1 * l + k;

            // load precomputed jacobian matrix terms
            c = jacobian[0][0][lindex] * jacobian[1][1][lindex] - jacobian[1][0][lindex] * jacobian[0][1][lindex];
            
            // load quadrature weight
            w = weights[0][k] * weights[1][l] * c;

            // add results to sumfactorization objects.
            // the terms in g are weighted by the quaderature weight.
            sumfact[0][lindex] = w;
            
        } // loop over 'k'
    } // loop over 'l'

} // end quadratureLoop


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // sizes
    mwSize m1;
    mwSize m2;
    
    // pointer objects
    double *sumfact[3];
    double *weights[2];
    double *jacobian[2][2];

    /* get the value of the scalar input  */
    m1 = mxGetScalar(prhs[0]);
    m2 = mxGetScalar(prhs[1]);
    
    /* create a pointer to the real data in the input matrix  */
    sumfact[0] = mxGetDoubles(mxGetCell(prhs[2], 0)); 
    sumfact[1] = mxGetDoubles(mxGetCell(prhs[2], 1)); 
    sumfact[2] = mxGetDoubles(mxGetCell(prhs[2], 2));

    weights[0] = mxGetDoubles(mxGetCell(prhs[3], 0)); 
    weights[1] = mxGetDoubles(mxGetCell(prhs[3], 1));

    jacobian[0][0] = mxGetDoubles(mxGetCell(prhs[4], 0)); 
    jacobian[1][0] = mxGetDoubles(mxGetCell(prhs[4], 1));
    jacobian[0][1] = mxGetDoubles(mxGetCell(prhs[4], 2)); 
    jacobian[1][1] = mxGetDoubles(mxGetCell(prhs[4], 3));

    /* call the computational routine */
    quadrature_loop(m1, m2,
        sumfact, weights,
                jacobian);
}
