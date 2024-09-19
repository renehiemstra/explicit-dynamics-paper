/*==========================================================
 * quadratureLoop.c - MATLAB External Interface
*
* Performs the quadrature loop to fill objects for use in 
* sum factorization, accompanying the provided matlab code.
*
* The calling syntax is:
*
*		void quadratureLoop(m1, m2, kappa, 
                            sumfact,
                            weights,
                            coordinates,
                            forcing,
                            displacement_ders,
                            jacobian)
*
* This is a MEX-file for MATLAB.
* Copyright 2007-2012 The MathWorks, Inc.
*
*========================================================*/

#include "mex.h"
#include "math.h"

double do_inverse_return_determinant(double A[2][2], double Ainv[2][2]){
    double c = 1.0 / (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
    Ainv[0][0] =  c * A[1][1]; Ainv[0][1] = -c * A[0][1];
    Ainv[1][0] = -c * A[1][0]; Ainv[1][1] =  c * A[0][0];
    return 1/c;
}

double norm(double a1, double a2){
    return sqrt(a1*a1 + a2*a2);
}

/* The computational routine */
void quadrature_loop(mwSize m1, mwSize m2, 
            double kappa,
            double *sumfact[3],
            double *weights[2],
            double *forcing,
            double *displacement_ders[2], 
            double *jacobian[2][2])
{
    mwSize k;
    mwSize l;
    mwSize lindex;

    double A[2][2], Ainv[2][2], DA[3][2];
    double Q[2][2], R[2][2];
    double x, y, w, c;
    double c_u, c_v;
    double c_grad[2], q[2], g[3];
    
    // loop over points in v-dir
    for (l=0; l<m2; l++) {
        // loop over points in u-dir
        for (k=0; k<m1; k++){
            // compute linear index
            lindex = m1 * l + k;

            // load precomputed jacobian matrix terms
            A[0][0] = jacobian[0][0][lindex]; A[0][1] = jacobian[0][1][lindex];
            A[1][0] = jacobian[1][0][lindex]; A[1][1] = jacobian[1][1][lindex];

            // compute inverse of jacobian matrix A and
            // return determinant of the mapping
            // c = qr_return_det(A, Q, R, Ainv);
            c = do_inverse_return_determinant(A, Ainv);
            
            // compute gradient of displacement
            q[0] = kappa * (Ainv[0][0] * displacement_ders[0][lindex] + Ainv[1][0] * displacement_ders[1][lindex]);
            q[1] = kappa * (Ainv[0][1] * displacement_ders[0][lindex] + Ainv[1][1] * displacement_ders[1][lindex]);

            // load quadrature weight
            w = weights[0][k] * weights[1][l] * c;

            // first two elements of 'g' relate to the regular gradient.
            // the third element of 'g' relates to taking the gradient of
            // the (1/c) term and the addition of the forcing term.
            g[0] = -Ainv[0][0] * q[0] - Ainv[0][1] * q[1];
            g[1] = -Ainv[1][0] * q[0] - Ainv[1][1] * q[1];
            g[2] = forcing[lindex];

            // add results to sumfactorization objects.
            // the terms in g are weighted by the quaderature weight.
            sumfact[0][lindex] = g[0] * w;
            sumfact[1][lindex] = g[1] * w;
            sumfact[2][lindex] = g[2] * w;
            
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
   
    // variables
    double kappa;
    
    // pointer objects
    double *sumfact[3];
    double *weights[2];
    double *forcing;
    double *displacement_ders[2];
    double *jacobian[2][2];

    /* get the value of the scalar input  */
    m1 = mxGetScalar(prhs[0]);
    m2 = mxGetScalar(prhs[1]);
    kappa = mxGetScalar(prhs[2]);
    
    /* create a pointer to the real data in the input matrix  */
    sumfact[0] = mxGetDoubles(mxGetCell(prhs[3], 0)); 
    sumfact[1] = mxGetDoubles(mxGetCell(prhs[3], 1)); 
    sumfact[2] = mxGetDoubles(mxGetCell(prhs[3], 2));

    weights[0] = mxGetDoubles(mxGetCell(prhs[4], 0)); 
    weights[1] = mxGetDoubles(mxGetCell(prhs[4], 1));

    forcing = mxGetDoubles(prhs[5]);

    displacement_ders[0] = mxGetDoubles(mxGetCell(prhs[6], 0)); 
    displacement_ders[1] = mxGetDoubles(mxGetCell(prhs[6], 1));

    jacobian[0][0] = mxGetDoubles(mxGetCell(prhs[7], 0)); 
    jacobian[1][0] = mxGetDoubles(mxGetCell(prhs[7], 1));
    jacobian[0][1] = mxGetDoubles(mxGetCell(prhs[7], 2)); 
    jacobian[1][1] = mxGetDoubles(mxGetCell(prhs[7], 3));

    /* call the computational routine */
    quadrature_loop(m1, m2, kappa,
        sumfact, weights,
            forcing, displacement_ders,
                jacobian);
}
