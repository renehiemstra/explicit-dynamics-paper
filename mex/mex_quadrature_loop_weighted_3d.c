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
                            jacobian,
                            hessian)
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


void set_jacobian(double A[3][3], double *jacobian[3][3], int lindex){
    for (mwSize j=0; j<3; j++) {
        for (mwSize i=0; i<3; i++) {
            A[i][j] = jacobian[i][j][lindex];
        }
    }
}

void set_hessian(double hessian[3][3][3], double *DA[6][3], int lindex){
    for (mwSize k=0; k<3; k++) {
        // main diagonal
        hessian[0][0][k] = DA[0][k][lindex];
        hessian[1][1][k] = DA[1][k][lindex];
        hessian[2][2][k] = DA[2][k][lindex];
        // above diagonal
        hessian[1][2][k] = DA[3][k][lindex];
        hessian[0][2][k] = DA[4][k][lindex];
        hessian[0][1][k] = DA[5][k][lindex];
        // below diagonal
        hessian[2][1][k] = DA[3][k][lindex];
        hessian[2][0][k] = DA[4][k][lindex];
        hessian[1][0][k] = DA[5][k][lindex];
    }
}

/* The computational routine */
void quadrature_loop(mwSize m1, mwSize m2, mwSize m3, 
            double kappa,
            double *sumfact[4],
            double *weights[3],
            double *forcing,
            double *displacement_ders[3], 
            double *jacobian[3][3],
            double *hessian[6][3])
{
    mwSize m;
    mwSize k;
    mwSize l;
    mwSize lindex;

    double A[3][3], Ainv[3][3], H[3][3][3];
    double w, c;
    double c_du[3], c_grad[3], q[3], g[4];

    // loop over points in v-dir
    for (m=0; m<m3; m++) {
        // loop over points in v-dir
        for (l=0; l<m2; l++) {
            // loop over points in u-dir
            for (k=0; k<m1; k++){
                // compute linear index
                lindex = m1 * m2 * m + m1 * l + k;
    
                // load precomputed jacobian matrix terms
                set_jacobian(A, jacobian, lindex);

                // load precomputed hessian terms
                set_hessian(H, hessian, lindex);

                // compute gradient of volume element (derivatives of jacobian)
                for (mwSize l=0; l<3; l++) {
                    c_du[l] = 0;
                    // first term in differentiating determinant
                    c_du[l] = c_du[l] +
                    H[l][0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) + 
                    A[0][0] * (H[l][1][1] * A[2][2] + A[1][1] * H[l][2][2]) -
                    A[0][0] * (H[l][2][1] * A[2][1] + A[1][2] * H[l][1][2]);
                    // second term in differentiating determinant
                    c_du[l] = c_du[l] -
                    H[l][1][0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + 
                    A[0][1] * (H[l][0][1] * A[2][2] + A[1][0] * H[l][2][2]) -
                    A[0][1] * (H[l][2][1] * A[2][0] + A[1][2] * H[l][0][2]);
                    // third term in differentiating determinant
                    c_du[l] = c_du[l] +
                    H[l][2][0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]) + 
                    A[0][2] * (H[l][0][1] * A[2][1] + A[1][0] * H[l][1][2]) -
                    A[0][2] * (H[l][1][1] * A[2][0] + A[1][1] * H[l][0][2]);
                }

                // compute inverse of jacobian matrix A and
                // return determinant of the mapping
                c = do_inverse_return_determinant(A, Ainv);

                // bring derivatives of jacobian to physical coordinates
                c_grad[0] = Ainv[0][0] * c_du[0] + Ainv[1][0] * c_du[1] + Ainv[2][0] * c_du[2];
                c_grad[1] = Ainv[0][1] * c_du[0] + Ainv[1][1] * c_du[1] + Ainv[2][1] * c_du[2];
                c_grad[2] = Ainv[0][2] * c_du[0] + Ainv[1][2] * c_du[1] + Ainv[2][2] * c_du[2];
                
                // compute gradient of displacement
                q[0] = kappa * (Ainv[0][0] * displacement_ders[0][lindex] + Ainv[1][0] * displacement_ders[1][lindex] + Ainv[2][0] * displacement_ders[2][lindex]);
                q[1] = kappa * (Ainv[0][1] * displacement_ders[0][lindex] + Ainv[1][1] * displacement_ders[1][lindex] + Ainv[2][1] * displacement_ders[2][lindex]);
                q[2] = kappa * (Ainv[0][2] * displacement_ders[0][lindex] + Ainv[1][2] * displacement_ders[1][lindex] + Ainv[2][2] * displacement_ders[2][lindex]);

                // load quadrature weight
                w = weights[0][k] * weights[1][l] * weights[2][m];
    
                // first two elements of 'g' relate to the regular gradient.
                // the third element of 'g' relates to taking the gradient of
                // the (1/c) term and the addition of the forcing term.
                g[0] = -Ainv[0][0] * q[0] - Ainv[0][1] * q[1] - Ainv[0][2] * q[2];
                g[1] = -Ainv[1][0] * q[0] - Ainv[1][1] * q[1] - Ainv[1][2] * q[2];
                g[2] = -Ainv[2][0] * q[0] - Ainv[2][1] * q[1] - Ainv[2][2] * q[2];
                g[3] = forcing[lindex] + ((q[0] * c_grad[0] + q[1] * c_grad[1] + q[2] * c_grad[2]) / c);
    
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
    double *hessian[6][3];

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

    hessian[0][0] = mxGetDoubles(mxGetCell(prhs[9], 0)); 
    hessian[1][0] = mxGetDoubles(mxGetCell(prhs[9], 1)); 
    hessian[2][0] = mxGetDoubles(mxGetCell(prhs[9], 2));
    hessian[3][0] = mxGetDoubles(mxGetCell(prhs[9], 3));
    hessian[4][0] = mxGetDoubles(mxGetCell(prhs[9], 4));
    hessian[5][0] = mxGetDoubles(mxGetCell(prhs[9], 5));

    hessian[0][1] = mxGetDoubles(mxGetCell(prhs[9], 6)); 
    hessian[1][1] = mxGetDoubles(mxGetCell(prhs[9], 7)); 
    hessian[2][1] = mxGetDoubles(mxGetCell(prhs[9], 8));
    hessian[3][1] = mxGetDoubles(mxGetCell(prhs[9], 9));
    hessian[4][1] = mxGetDoubles(mxGetCell(prhs[9], 10));
    hessian[5][1] = mxGetDoubles(mxGetCell(prhs[9], 11));

    hessian[0][2] = mxGetDoubles(mxGetCell(prhs[9], 12)); 
    hessian[1][2] = mxGetDoubles(mxGetCell(prhs[9], 13)); 
    hessian[2][2] = mxGetDoubles(mxGetCell(prhs[9], 14));
    hessian[3][2] = mxGetDoubles(mxGetCell(prhs[9], 15));
    hessian[4][2] = mxGetDoubles(mxGetCell(prhs[9], 16));
    hessian[5][2] = mxGetDoubles(mxGetCell(prhs[9], 17));

    /* call the computational routine */
    quadrature_loop(m1, m2, m3, 
        kappa,
        sumfact, 
        weights,
        forcing, 
        displacement_ders,
        jacobian, 
        hessian);
}
