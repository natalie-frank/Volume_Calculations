/*
maps the 90 degree torus to the 60 degree torus
for measuring distances on the 60 degree torus using represenations on the 90
degree torus

The calling syntax is:
 
 		t = tor_90_to_tor_60(z)

where 
       z = vector of particle positions on the 90 degree torus (2n x 1 column vector)
       t = vector of particle positions on the 60 degree torus (2n x 1 column vector)
*/


#include "mex.h"
#include "math.h"




void torus_conversion(double* z, int n, double* t)
{
    int i;
    double zm[2*n];
    double val;
    for(i=0;i<2*n;i++)
    {
        val=z[i];
        while(val<0|| val>=1)
        {
            if(val<0)
            {
                val++;
            }
            else
            {
                val--;
            }
        }
        zm[i]=val;
        t[i]=val;
    }
    for( i=0;i< 2*n; i++)
    {
        if(i % 2==0)
        {

            t[i]=zm[i]+zm[i+1]/2;
        }
        else
        {
            t[i]= sqrt(3)/2* zm[i];
        }
    }
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *z;               /* the input vector on the 90 degree torus */
    double *t;               /* the output vector on the 60 degree torus */
    
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
    z = mxGetPr(prhs[0]);  /* pointer to z */
    n = mxGetM(prhs[0])/2;   /* number of particles */
    

    /* Create the output matrix, initialized to 0 */
    plhs[0] = mxCreateDoubleMatrix((mwSize)2*n,1,mxREAL);
    
    /* get a pointer to the real data in the output matrix */
    t = mxGetPr(plhs[0]);

    /* call the computational routine */
    torus_conversion(z, n, t);
}