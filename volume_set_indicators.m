
        %%inputs%%
%N is the number of samples

%burn_in is the burn_in time

%up: up is true if we compute the volume of our desired shape by comparing
%to a larger shape and false if we compar to a smaller shape

%increasing: increasing is true if our shape is given by the set of points
%which satisfy g_i(x)\leq k_i and false the shape is given by g_i(x)\geq
%k_i. Note that increasing is true if the set of shapes given by M(x)\leq
%increases with m and is false if it decreases with m

%x0 is the initial point the method

%next: This is a mechanism that allows varying the coordinates in a gibbs
%sampler. Specifically, next is a function that takes the coordinates for
%which the method just make a propsal as input and outputs the next
%coordinates for which the MCMC method should make a proposal

%start_coordinates: start_coordinates is the first set of coordinates we
%use in our Gibbs sampler

%proposal_function: this function takes in two arguments (x,k), the current
%point and the coordinates from which to sample. it outputs the proposal
%move in the MCMC method

%H: H takes in 3 arguments (y,x,k). y is the proposal move, x is the
%current point, and k is the coordinates that were sampled to get from k to
%y. H returns the P(y->x)/P(x->y).

%m: The shape(s) of which we want to find the volume are given by M(x)\leq
%m if increasing is true, and M(x)\geq m if increasing is false. m can be a
%list

%M takes in a vector x outputs the 'smallest' shape which contains x. Formally,
%M outputs \max_i g_i(x)/k_i if increasing is true, and \min_i g_i(x)/k_i if increasing is false 
%reference_volume is the volume of the largest shape 

%return_samples: if return_samples is true, we return the samples we drew along with the 
%random variables in the numerator/denominator of our sum, and which samples
%were accepted. Otherwise, we return []

%W_importance is a function that takes in a sample x and outputs the
%importance weight at x

%W_estimator weights our samples when calculating the ratio of our two volumes.
%Explicitly, in the 'volume_ratio' function, the volume ratio is esimated 
%as \sum_{i=1}^N indicator_numerator(x,m)* W_estimator(x) / \sum_{i=1}^N
%indicator_denominator(x,m) W_estimator(x)
%('indicator_numerator' and 'indicator_denominator' are two other inputs to
%volume_ratio)


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator

%accepted is whether or not the proposal was accepted at each monte carlo step
function [ratio,samples,numerator_rv,denominator_rv,accepted]  = volume_set_indicators(N,burn_in,up,increasing,x0, next,start_coordinates, proposal_function,H,m,M, return_samples,W_importance,W_estimator)
    len=length(m);    
    if increasing
        increasing_multiplier=1;
    else
        increasing_multiplier=-1;
    end
    if up    
        numerator_indicator=@(x,mm)increasing_multiplier*M(x)<=increasing_multiplier*mm;
        denominator_indicator=@(x,mm) ones(len,1); 
    else
        numerator_indicator=@(x,mm) ones(len,1); 
        denominator_indicator=@(x,mm)increasing_multiplier*M(x)<=increasing_multiplier*mm;
    end
    
    [ratio,samples,numerator_rv,denominator_rv,accepted]=volume_ratio(N,burn_in,x0, next,start_coordinates,proposal_function,H,m,W_importance,W_estimator, numerator_indicator,denominator_indicator,return_samples); 
    
end

