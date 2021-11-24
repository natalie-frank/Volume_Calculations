/*==========================================================
 * distx2d.c -- calculates distance matrix between points, 
 *              assuming periodic boundary conditions with 
 *              side length 1. 
 *
 * The calling syntax is:
 *
 *		R = distx2dP(x)
 *
 * where 
 *      x = vector of particle positions (2n x 1 column vector)
 *      R (output) = n x n distance matrix
 *
 *
 * Note that double arrays are indexed as (row) + (Nrows)*(col)
 *
 *========================================================*/
/* $Revision: 1.1.10.4 $ */

#include "mex.h"
#include "math.h"

/* The computational routine */

void distx2d(double *x, int n, double *R){
    
    int i,j,k;  /* for loops */
    double d3;  /* holds distance */
    double dx, dy;  /* holds distances in x, y */
	
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++) {
            
            dx = x[2*i]-x[2*j];
            dy = x[2*i+1]-x[2*j+1];
            
            if(dx > 0.5) dx = dx - 1;
            if(dx < -0.5)  dx = dx + 1;
            if(dy > 0.5) dy = dy - 1;
            if(dy < -0.5)  dy = dy + 1;
            
            d3 = sqrt(dx*dx + dy*dy);      

            R[i + n*j] = d3;
            R[j + n*i] = d3;
		}
	}
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *x;               /* vector x */
    double *R;               /* output matrix */
    
    int n;                   /* number of particles */
    
    size_t mrows,ncols;      /* for error checking */

    /* Check inputs */
    
    /* Check for proper number of arguments. */
    if(nrhs!=1) {
      mexErrMsgIdAndTxt( "MATLAB:distx:invalidNumInputs",
              "One input required.");
    } else if(nlhs>1) {
      mexErrMsgIdAndTxt( "MATLAB:distx:maxlhs",
              "Too many output arguments.");
    }
    
    /* The input must be a noncomplex vector double.*/
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
        !(ncols==1) ) {
        mexErrMsgIdAndTxt( "MATLAB:distx:inputNotRealColumnVectorDouble",
              "Input must be a noncomplex column vector double.");
    }

    /* Number of rows must be a multiple of 2 */
    if( !(mrows % 2 == 0) ) {
        mexErrMsgIdAndTxt( "MATLAB:distx:inputNotParticlePositions",
              "Vector length must be divisble by 2.");
    }
    
    /* Get input variables */
    x = mxGetPr(prhs[0]);  /* pointer to x */
    n = mxGetM(prhs[0])/2;   /* number of particles */
    

    /* Create the output matrix, initialized to 0 */
    plhs[0] = mxCreateDoubleMatrix((mwSize)n,(mwSize)n,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    R = mxGetPr(plhs[0]);

    /* call the computational routine */
    distx2d(x, n, R);

}
