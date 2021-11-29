/*

x is the disk locations, odd indices the x-axis, even indices the y-axis
k is the index of the disk we want to move. k is passed as a double and then truncated
r is the radius of the disks
M is a PSD matrix defining the inner product on our domain
finds the maximum distance disk k could travel before hitting another
disk, assuming it moves off the edge of the torus at most three times

calling syntax:

max_dist=maximum_move(x,k,r, theta,A)


*/


#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "float.h"
#include "limits.h"





//computes the product A^TA for a 2x2 matrix A and stores it in M
//M must be a a double array length 4
//Assumes that A,M are in the order A[0]=A_{1,1} A[1]=A_{2,1} A[2]=A_{1,2} A[3]=A_{2,2}
void get_M(double* A, double *M)
{
    M[0]=A[0]*A[0]+A[1]*A[1];
    M[1]=A[0]*A[2]+A[1]*A[3];
    M[2]=A[0]*A[2]+A[1]*A[3];
    M[3]=A[2]*A[2]+A[3]*A[3];
}


//returns the minimum of two doubles
double min(double a, double b)
{
    if(a<=b)
    {
        return a;
    }
    else
    {
        return b;
    }
}

//returns the maximum of two doubles
double max(double a, double b)
{
    return -1*min(-1*a,-1*b);
}

//returns the 2-norm of an array with 2 elements
double norm(double *x)
{
    double out;// the output
    out=sqrt(x[0]*x[0]+x[1]*x[1]);
    return out;
}


//If x,y are arrays length 2 and M is a 2x2 matrix, then this returns x^T M y
//(Or maybe y^T M x, but in our case it shouldn't matter as M should be symmetric)
double mtrx_dot(double *x,double *M, double *y)
{
    double *ptr=&M[0];
    double out;//the output
    int i=0;
    int j=0;
    int dim=2;
    for(i=0;i<dim;i++)
    {
        for(j=0;j<dim;j++)
        {
            out=out+y[i]*x[j]* *ptr;
            if(i<dim-1|| j<dim-1) ptr=ptr+1;
        }
    }
    
    return out;
    
}

//returns the dot product of x and y for 2x1 vectors
double dot(double *x ,double *y)
{
    return x[0]*y[0]+x[1]*y[1];
}





//Finds the (t,s) for which the lines x+tv and y+sw intersect, assuming v,w are not parallel
void line_intersection(double*x,double*v,double*y,double*w,double*intersection)
{
    double diff[2];
    diff[0]=x[0]-y[0];
    diff[1]=x[1]-y[1];
    double det=w[1]*v[0]-w[0]*v[1];
    if(det==0)//the lines don't intersect or they are equal. If they are equal the maximum distance is 3-2*horizontal_diff, and we already set L to be less than or equal to this quantity.
    {
        //we set the solutions to the minimum (negative) integer, as this won't affect the program later on
        intersection[0]=INT_MIN;
        intersection[1]=INT_MIN;
    }
    /*else solution is given by
     *x-y=[ w(0) v(0)  * [s
     *      w(1) v(1)]    t]
     */

    else
    {
        double T_inv[2][2]={ { -1*v[1], v[0]}, //if we call the matrix in the comment above T, T_inv is its inverse
                         {-1*w[1], w[0]}};

        double s_star=(T_inv[0][0]*diff[0]+T_inv[0][1]*diff[1])/det;
        double t_star=(T_inv[1][0]*diff[0]+T_inv[1][1]*diff[1])/det;

        intersection[0]=t_star;
        intersection[1]=s_star;
    }

}



void max_dist(double* x, int* k, double r, double theta, double* A, double* output, int n)
//void max_dist(double* x, int k, double r, double theta, double* A, bool* output, int n)
{
    double e1[2]={1,0};
    double e2[2]={0,1};
    double L;
    //double locs[9*(2*n-2)];//list wil ocntain all translations of the disks we need to consider
    double x0[2];//the current disk
    //double xx[2*n-2];//an array of all the other disks
    //int m=9*(2*n-2);//the number of disks under consideration
    double d[2]={cos(theta),sin(theta)};//the direction of motion
    double M[4];//for A^TA
    double disk[2];// the current disk under consideration
    double diff[2];// for storing the difference between disk positions
    double dist; //for storing the distance between 2 disks
    double a;//stores the constant a in our quadratic ax^2+bx+c
    double b_prime;//stores -b/2 for our quadratic ax^2+bx+c
    double discriminant_div;//for storing the discriminant divided by 4a
    double t; // for storing the point on x_0+td where the line gets within 2r of a disk center
    int i;//index variable
    int j;//index variable
    int ell;//index variable
    double v[2];// for storing vectors during the computation of L
    double pt[2];//for storing the four vertices of our parallelogram
    
    //the current disk
    x0[0]=x[k[0]];
    x0[1]=x[k[1]];

    

    //finding M and the inverse of A
    get_M(A,M);
    
    
    
   


    //The largest possible distance occurs on the diagonals of the parallelogram. We compute the magnitude of the two diagonals and take the maximum.
    //We multiply by 3 because we consider a 3x3 box of parallelograms
    v[0]=A[0]+A[2];
    v[1]=A[1]+A[3];
    L=3*norm(v);
    v[0]=A[0]-A[2];
    v[1]=A[1]-A[3];
    L=max(L,3*norm(v));
    

    
    double col1[2]={A[0], A[1]};
    double col2[2]={A[2], A[3]};
    double cos_alpha=dot(col1,col2)/(norm(col1)*norm(col2));
    double sin_alpha=sqrt(1-cos_alpha*cos_alpha);
    double r_mapped=r/sin_alpha;
    
    double t_max=3-4*r_mapped;//3 for 3 copies of the torus;


    double xs [4][2]={
                        {2-2*r_mapped,2-2*r_mapped},
                        {-1+2*r_mapped,2-2*r_mapped},
                        {-1+2*r_mapped,-1+2*r_mapped},
                        {2-2*r_mapped,-1+2*r_mapped}
                    };
    double vs [4][2]={    //moving around the corners of the parallelogram in a counterclockwise direction
                        {-1,0},
                        {0,-1},
                        {1,0},
                        {0,1}
                    };
    double intersection[2];
    for(i=0;i<=3; i++)
    {   v[0]=vs[i][0];
        v[1]=vs[i][1];
        pt[0]=xs[i][0];
        pt[1]=xs[i][1];
        line_intersection(x0,d,pt,v,intersection);
        if(intersection[0]>=0&& intersection[1]<=t_max&&intersection[1]>=0)
        {
            L=min(L,intersection[0]);
            
        }

        

    }




    *output=L;
    
    double discriminant;
    double c;
    
    for(i=-1; i<=1; i++)
    {
        for(j=-1; j<=1; j++)
        {
            for(ell=0; ell<n; ell++)
            {
                if(ell!=k[1]/2)
                {
                    disk[0]=x[2*ell]+i;
                    disk[1]=x[2*ell+1]+j;
                    diff[0]=disk[0]-x0[0];
                    diff[1]=disk[1]-x0[1];
                    b_prime=mtrx_dot(diff,M,d);// this is -b/2

                    
                    if(b_prime>0) //if mtrx_dot(diff,M,d)<=0, the closest point between y and x+td is for negative (or 0) t. Thus the distance only increases if t>0, x0 will not run into y when moving in direction d
                    {
                        a=mtrx_dot(d,M,d);
                        c=mtrx_dot(diff,M,diff)-4*r*r;

                        if(4*r*r>mtrx_dot(diff,M,diff) || b_prime*b_prime>=a*c)//we compare before subtracting (twice) due to conditioning
                        {
                            discriminant=b_prime*b_prime-a*c;
                             //the if contition above implies that the discriminant should be non-negative. However, the discriminant may end up below zero due to rounding error. If this occurs, we round up to 0.
                            discriminant=max(discriminant,0);

                            t=(b_prime-sqrt(discriminant))/a;

                            t=max(t,0);//we add this statement to guard agains numerical error-- we should always have t>=0 as b>= discriminant in our situation
                            *output=min(t,*output);
                        }
                    }
                }
            }
        }
    }
}










void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double* x;//the disk locations
    int k[2];//the indices of the disk we want to move
    double* r;//the radius of the disks
    double* theta;// the angle by which we move the disk
    double* A;// the PSD matrix defining the inner product on our domain
  
    double* mx_dist;//the output

    int n;//the number of disks

    int mrows;//for keeping track of the number of rows in inputs
    int ncols;//for keeping track of the number of columns in inputs

       /* Check for proper number of arguments. */
    if(nrhs!=5) 
    {
      mexErrMsgIdAndTxt( "MATLAB:max_dist:invalidNumInputs",
              "Five inputs required.");
    } 




    //Get and check the first input
        /* The input must be a noncomplex vector double.*/
    mrows = mxGetM(prhs[0]);
    ncols = mxGetN(prhs[0]);
    if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0])) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotRealMatrixDouble",
              "First Input must be a noncomplex matrix double.");
    }

    /* Number of rows must be a multiple of 2 */
    if( !(mrows % 2 == 0) ) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotParticlePositions",
              "Number of matrix rows must be divisble by 2.");
    }
    
    /* Get the first input variable */
    x = mxGetPr(prhs[0]);  /* pointer to x */
    n = mrows/2;   /* number of particles */




//Get and check the second input
    if(!mxIsDouble(prhs[1]))
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:diskIndicesNotInteger", "The disk inices must be two Integers");
    }

    
    mrows = mxGetM(prhs[1]);
    ncols = mxGetN(prhs[1]);

    if(!(mrows==1 && ncols== 2)&& !(mrows==2 && ncols==1)) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotIntegerVectorLengthTwo",
              "Second input must be a integer vector length 2.");
    }
    double* pointer=mxGetPr(prhs[1]);
    k[0]=(int)pointer[0];
    k[1]=(int)pointer[1];
    if(k[0]<1||k[1]>2*n|| k[1]!=k[0]+1||k[1]%2!=0)
    {         
                mexErrMsgIdAndTxt( "MATLAB:max_dist:inputOutOfRange",
              "Second input must be two consecutive integers between zero and twice the number of disks , the first of which is odd.");
    }
    k[0]=k[0]-1;//convert 1-based indexing to 0-based indexing
    k[1]=k[1]-1;


    //Get and check the third input

    if(!mxIsDouble(prhs[2]))
    {
                mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotRealDouble",
              "Third Input must be a noncomplex double.");
    }

    r=mxGetPr(prhs[2]);



    mrows = mxGetM(prhs[2]);
    ncols = mxGetN(prhs[2]);
    if(mrows!=1 || ncols!= 1) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotScalarDouble",
              "The third input must be a scalar.");
    }




    //get and check the fourth input

    if(!mxIsDouble(prhs[3]))
    {
                mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotRealDouble",
              "Fourth Input must be a noncomplex double.");
    }

    theta=mxGetPr(prhs[3]);



    mrows = mxGetM(prhs[3]);
    ncols = mxGetN(prhs[3]);
    if(mrows!=1 || ncols!= 1) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotScalarDouble",
              "The fourth input must be a scalar.");
    }






       //Get and check the fifth input
        /* The input must be a noncomplex vector double.*/
    mrows = mxGetM(prhs[4]);
    ncols = mxGetN(prhs[4]);
    if( !mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) 
    {
        mexErrMsgIdAndTxt( "MATLAB:max_dist:inputNotRealMatrixDouble",
              "Fifth Input must be a noncomplex matrix double.");
    }
    if(mrows!=2|| ncols!=2)
    {
                mexErrMsgIdAndTxt( "MATLAB:max_dist:incorrectMatrixDimensions",
              "Fifth Input must be 2x2 matrix.");
    }
    A=mxGetPr(prhs[4]);



    //Create the output

    plhs[0]=mxCreateDoubleScalar(0);
    mx_dist=mxGetPr(plhs[0]);

    //plhs[0]=mxCreateLogicalMatrix(3,3);

    //bool* output=mxGetLogicals(plhs[0]);



    max_dist(x, k, *r, *theta, A, mx_dist, n);
    //max_dist(x, k, *r, *theta, A, output, n);


}
