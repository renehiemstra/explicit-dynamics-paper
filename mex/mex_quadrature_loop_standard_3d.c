/*==========================================================
 * quadratureLoop.c - MATLAB External Interface
*
* Performs the quadrature loop to fill objects for use in 
* sum factorization, accompanying the provided matlab code.
*
* The calling syntax is:
*
*		void quadratureLoop(m1, m2, m3, kappa, 
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

void  printmat(double A[3][3]){
    printf("[%0.3f   %0.3f   %0.3f\n",  A[0][0], A[0][1], A[0][2]);
    printf(" %0.3f   %0.3f   %0.3f\n",  A[1][0], A[1][1], A[1][2]);
    printf(" %0.3f   %0.3f   %0.3f]\n", A[2][0], A[2][1], A[2][2]);
}

double det(double a, double b, double c, double d){
    return a * d - b * c;
}

double do_inverse_return_determinant(double A[3][3], double Ainv[3][3]){
    // compute transpose of the cofactor matrix
    for (mwSize i=0; i<3; i++) {
        for (mwSize j=0; j<3; j++) {
            Ainv[j][i] = det(A[(i+1)%3][(j+1)%3], A[(i+1)%3][(j+2)%3], A[(i+2)%3][(j+1)%3], A[(i+2)%3][(j+2)%3]);
        }
    }
    // compute determinant of the matrix
    double d = 0.0;
    for (mwSize i=0; i<3; i++) {
        d = d + A[i][0] * Ainv[0][i];
    }
    // mutiply transpose cofactor matrix by inverse of determinant
    // to get the inverse matrix
    for (mwSize j=0; j<3; j++) {
        for (mwSize i=0; i<3; i++) {
            Ainv[i][j] = Ainv[i][j] / d;
        }
    }
    return d;
}

double norm(double a1, double a2, double a3){
    return sqrt(a1*a1 + a2*a2 + a3*a3);
}

/* The computational routine */
void quadrature_loop(mwSize m1, mwSize m2, mwSize m3, 
            double kappa,
            double *sumfact[4],
            double *weights[3],
            double *forcing,
            double *displacement_ders[3], 
            double *jacobian[3][3])
{
    mwSize m;
    mwSize k;
    mwSize l;
    mwSize lindex;

    double A[3][3], Ainv[3][3];
    double w, c;
    double q[3], g[4];

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

                // compute inverse of jacobian matrix A and
                // return determinant of the mapping
                c = do_inverse_return_determinant(A, Ainv);
                
                // compute gradient of displacement
                q[0] = kappa * (Ainv[0][0] * displacement_ders[0][lindex] + Ainv[1][0] * displacement_ders[1][lindex] + Ainv[2][0] * displacement_ders[2][lindex]);
                q[1] = kappa * (Ainv[0][1] * displacement_ders[0][lindex] + Ainv[1][1] * displacement_ders[1][lindex] + Ainv[2][1] * displacement_ders[2][lindex]);
                q[2] = kappa * (Ainv[0][2] * displacement_ders[0][lindex] + Ainv[1][2] * displacement_ders[1][lindex] + Ainv[2][2] * displacement_ders[2][lindex]);

                // load quadrature weight
                w = weights[0][k] * weights[1][l] * weights[2][m] * c;
    
                // first two elements of 'g' relate to the regular gradient.
                // the third element of 'g' relates to taking the gradient of
                // the (1/c) term and the addition of the forcing term.
                g[0] = -Ainv[0][0] * q[0] - Ainv[0][1] * q[1] - Ainv[0][2] * q[2];
                g[1] = -Ainv[1][0] * q[0] - Ainv[1][1] * q[1] - Ainv[1][2] * q[2];
                g[2] = -Ainv[2][0] * q[0] - Ainv[2][1] * q[1] - Ainv[2][2] * q[2];
                g[3] = forcing[lindex];
    
                // add results to sumfactorization objects.
                // the terms in g are weighted by the quaderature weight.
                sumfact[0][lindex] = g[0] * w;
                sumfact[1][lindex] = g[1] * w;
                sumfact[2][lindex] = g[2] * w;
                sumfact[3][lindex] = g[3] * w;
                
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
   
    // variables
    double kappa;
    
    // pointer objects
    double *sumfact[4];
    double *weights[3];
    double *forcing;
    double *displacement_ders[3];
    double *jacobian[3][3];

    /* get the value of the scalar input  */
    m1 = mxGetScalar(prhs[0]);
    m2 = mxGetScalar(prhs[1]);
    m3 = mxGetScalar(prhs[2]);
    kappa = mxGetScalar(prhs[3]);
    
    /* create a pointer to the real data in the input matrix  */
    sumfact[0] = mxGetDoubles(mxGetCell(prhs[4], 0)); 
    sumfact[1] = mxGetDoubles(mxGetCell(prhs[4], 1)); 
    sumfact[2] = mxGetDoubles(mxGetCell(prhs[4], 2));
    sumfact[3] = mxGetDoubles(mxGetCell(prhs[4], 3));

    weights[0] = mxGetDoubles(mxGetCell(prhs[5], 0)); 
    weights[1] = mxGetDoubles(mxGetCell(prhs[5], 1));
    weights[2] = mxGetDoubles(mxGetCell(prhs[5], 2));

    forcing = mxGetDoubles(prhs[6]);

    displacement_ders[0] = mxGetDoubles(mxGetCell(prhs[7], 0)); 
    displacement_ders[1] = mxGetDoubles(mxGetCell(prhs[7], 1));
    displacement_ders[2] = mxGetDoubles(mxGetCell(prhs[7], 2));

    jacobian[0][0] = mxGetDoubles(mxGetCell(prhs[8], 0)); 
    jacobian[1][0] = mxGetDoubles(mxGetCell(prhs[8], 1));
    jacobian[2][0] = mxGetDoubles(mxGetCell(prhs[8], 2));
    jacobian[0][1] = mxGetDoubles(mxGetCell(prhs[8], 3)); 
    jacobian[1][1] = mxGetDoubles(mxGetCell(prhs[8], 4));
    jacobian[2][1] = mxGetDoubles(mxGetCell(prhs[8], 5));
    jacobian[0][2] = mxGetDoubles(mxGetCell(prhs[8], 6)); 
    jacobian[1][2] = mxGetDoubles(mxGetCell(prhs[8], 7));
    jacobian[2][2] = mxGetDoubles(mxGetCell(prhs[8], 8));

    /* call the computational routine */
    quadrature_loop(m1, m2, m3, 
        kappa,
            sumfact, weights,
                forcing, displacement_ders,
                    jacobian);
}
