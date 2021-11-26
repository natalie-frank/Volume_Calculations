/*==========================================================
 * dist_pairs_60_on_90.c -- calculates distance matrix between points on 60 degree torus using the representation of the 60 degree torus on the 90 the degree torus (ie, the unit square). 
 *              We assume that this 60 degree torus has periodic boundary conditions with 
 *              side length 1. 
 *
 * The calling syntax is:
 *
 *		R = dist_pairs_60_on_90(x)
 *
 * where 
 *      x = vector of particle positions (2n x 1 column vector) odd entries are the x-coordinates and even entries are the y-coordinates
 *      R (output) = n x n distance matrix
 *
 *
 * Note that double arrays within C are indexed as (row) + (Nrows)*(col)
 *
 *
 *========================================================*/


#include "mex.h"
#include "math.h"

/* The computational routines */
void torus_conversion(double* z, int n, double* t)
{
    int i;
    double zm[2*n];// will store z mod 1
    double val;
    //finding z mod 1
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
    }
    //converting to the 60 degree torus
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


void distx2d(double *x, int n, double *R){
    
    int i,j,k;  /* for loops */
    double dist;  /* distance */
    double d1, d2, d3, d4;  /* distances for different translations */
    double rx, ry;  /* separation */
    double l1x,l1y,l2x,l2y;  /* lattice vectors */
    double s1,s2;  /* signs */
    
    l1x = 1;
    l1y = 0;
    l2x = 0.5;
    l2y = sqrt(3)/2;
	
	for(i=0; i<n; i++){
		for(j=i+1; j<n; j++) {
            
            rx = x[2*i]-x[2*j];
            ry = x[2*i+1]-x[2*j+1];
 
            s1 = 0;
            s2 = 0;
            
            if(rx*l1x+ry*l1y > 0)   s1 = -1;
            if(rx*l1x+ry*l1y < 0)  s1 = 1;
            if(rx*l2x+ry*l2y > 0)   s2 = -1;
            if(rx*l2x+ry*l2y < 0)  s2 = 1;
            /* changed 0.5 to 0 in the above, because of problems
             * with the diagonal term */ 

            dist = rx*rx + ry*ry;
            
            /*
            printf("rx = %f, ry = %f\n",rx, ry);
            printf("r.l1 = %f\n",rx*l1x+ry*l1y);
            printf("r.l2 = %f\n",rx*l2x+ry*l2y);
            printf("dist = %f\n",dist);
            printf("s1 = %f, s2 = %f\n",s1,s2);
             */
            
            if(fabs(s1) > 0) {
                d2 = (rx+s1*l1x)*(rx+s1*l1x) + (ry+s1*l1y)*(ry+s1*l1y);
                dist = fmin(dist,d2);
                
                //printf("d2 = %f\n",d2);
            } 
            if(fabs(s2) > 0) {
                d3 = (rx+s2*l2x)*(rx+s2*l2x) + (ry+s2*l2y)*(ry+s2*l2y);
                dist = fmin(dist,d3);
                
                //printf("d3 = %f\n",d3);
            } 
            if(fabs(s1) > 0 && fabs(s2) > 0) {
                d4 = (rx+s1*l1x+s2*l2x)*(rx+s1*l1x+s2*l2x) + (ry+s1*l1y+s2*l2y)*(ry+s1*l1y+s2*l2y);
                dist = fmin(dist,d4);
                
                //printf("d4 = %f\n",d4);
            } 

            dist = sqrt(dist);

            R[i + n*j] = dist;
            R[j + n*i] = dist;
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

    /* call the computational routines */

    double t[2*n];
    torus_conversion(x, n, t);
    distx2d(t, n, R);

}
