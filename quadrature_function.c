/*
nodes is nodes for quadrature
IMPORTANT ASSUMPTION: nodes is sorted

w is weights for the quadrature
the function uses the trapezoidal rule to numerically integrate w from
min(nodes) to y
weights are the quadrature weights
trapezoidal quadrature is assumed for the two points with nodes(m)\leq
y\leq nodes(m+1)

assumes that y>=min(nodes)

calling syntax:

y=quadrature_function(nodes,weights,w,y)
*/

#include "mex.h"



//Find the minimum of 2 integers

int min(int a1, int a2)
{
    if(a1<=a2)
    {
        return a1;
    }
    else
    {
        return a2;
    }
    
}

//Find the maximum of two integers
int max(int a1, int a2)
{
    if(a1>=a2)
    {
        return a1;
    }
    else
    {
        return a2;
    }
    
}


//The computational routine

void quad_f(double* nodes, double* weights, double* w, double* y, int len, double* integral)
{
    
    int m=0;//the integer m for which nodes[m]<=y<=nodes[m+1]
    double weight_m;// newly computed mth weight
    double weight_mp1;// newly computed m+1th weight
    double delta;//difference between nodes[m+1] and nodes[m]
    double delta_prev;//difference between nodes[m] and nodes[m-1]

    //Finding the integer m for which nodes[m]<=y<nodes[m+1]
    if(nodes[len-1]<=*y)// at the last node, this isn't possible so we take m=len-1;
    {
        m=len-1;
    }
    else
    {
        while(*y>=nodes[m+1])//note that this condition implies that nodes[m+1]>nodes[m], so lower down delta!= 0
        {
            m++;
        }
        
    }
    


    //if we're at the max m, then the last 2 weights don't need to be changed, but we still sum from 0 to m+1, so we decrement m by 1
    if(m==len-1)
    {
        m--;
        weight_m=weights[m];
        weight_mp1=weights[m+1];
    }
    //otherwise, we re-compute the last two weights
    else
    {
        delta= nodes[m+1]-nodes[m];
        if(m>0)
        {
            delta_prev=nodes[m]-nodes[m-1];
            weight_m=(*y-nodes[m-1])/2;
        }
        else
        {
            weight_m=(*y-nodes[0])/2;
            delta_prev=0;
        }
        double p=(nodes[m+1]- *y)/delta;
        
        weight_mp1=delta/2*(1-p);
    }
    
    //now we take the dot product of weights and the function values w to get the value of the integral

    // up to m-1, we use the passed in weights
    int i;
    for(i=0; i<=m-1;i++)
    {
        *integral+= weights[i] * w[i];
    }
    // for m and m+1, we use the newly computed weights
    *integral+= weight_m * w[m];
    *integral+= weight_mp1 * w[m+1];
}

//the mex gateway function
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *nodes;   // the quadrature nodes
    double *weights; // the quadrature weights
    double *w; // the function evaluated at the nodex
    double *y; // the point we want to integrate to
    double *integral; // the resulting integral


    // The sizes of our input variables
    /*size_t nodes_rows, nodes_cols;
    size_t weights_rows, weights_cols;
    size_t w_rows, w_cols;
    size_t y_rows, y_cols;*/

    // the number of rows and columns in each of our input variables
    size_t rows[4];
    size_t cols[4];

    //vector lengths for the first three inputs
    size_t vect_length[3];

    int i;//loop variable

    //checking the inputs
     /* Check for proper number of arguments. */
    if(nrhs!=4) 
    {
        mexErrMsgIdAndTxt( "MATLAB:distx:invalidNumInputs",
              "Four inputs required.");
    }

    /* All inputs must be noncomplex vector doubles.*/
    for(i=0; i<4; i++ )
    {
        if( !mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]))
        {
            mexErrMsgIdAndTxt( "MATLAB:distx:inputNotRealColumnVectorDouble",
              "All inputs must be noncomplex doubles.");
        } 
    }
    for(i=0; i<4; i++)
    {
        rows[i]=mxGetM(prhs[i]);
        cols[i]=mxGetN(prhs[i]);
    }

    //check that the last input is a 1x1 matrix
    if( rows[3]!=1 || cols[3]!=1)
    {
        mexErrMsgIdAndTxt("MATLAB:quad_func:lastInputNotScalar", "Last input must be a 1x1 double");
    }

    // check that each of first three inputs is a 1xn or nx1 vector, and then puts the vector length in vect_length
    for(i=0;i<3; i++)
    {
        if(rows[i]!=1 && cols[i]!=1)
        {
            mexErrMsgIdAndTxt("MATLAB:quad_func:inputDimensionsIncorrect","The first three inputs must be nx1 or 1xn vectors");
        }
        vect_length[i]=max(rows[i],cols[i]);
    }
    
    //Check that the first three inputs have the same length
    int len=vect_length[0];
    if(vect_length[1]!=len || vect_length[2]!=len)
    {
        mexErrMsgIdAndTxt("MATLAB:quad_func:inputDimensionsIncorrect","The first three inputs must be vectors of the same length");
    }
    //Get the input variables
    nodes = mxGetPr(prhs[0]);
    weights = mxGetPr(prhs[1]);
    w = mxGetPr(prhs[2]);
    y = mxGetPr(prhs[3]);

    //Check that y is in the range of nodes
    if(*y<nodes[0])
    {
        mexErrMsgIdAndTxt("MATLAB:quad_fund:quadraturePointOutOfRange","The quadrature point y is less than the minimum node");
    }
    /*else if(*y>nodes[len-1])
    {
        mexErrMsgIdAndTxt("MATLAB:quad_fund:quadraturePointOutOfRange","The quadrature point y is greater than the maximum node");
    }*/

    //Get a pointer to the output
    plhs[0]=mxCreateDoubleMatrix( 1,1, mxREAL);
    integral=mxGetPr(plhs[0]);

    //call the computational routine
    quad_f(nodes, weights, w, y, len, integral);
}