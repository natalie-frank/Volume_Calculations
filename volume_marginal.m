
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
%point and the coordinates from which to sample.
%There are two options for the output:
%1. If using MCMC to sample W(M(x)), it outputs the proposal move in the method. 
%H (the next parameter) must be a funtion which outputs P(y\to x)/P(x\to y)
%2. it can directly sample from the distribution proportional to W(M(x)). H must be []

%H: There are 2 options for H
%1. If using MCMC to sample from W(M(x)), H takes in 3 arguments (y,x,k). y is the proposal 
%move, x is the current point, and k is the coordinates that were sampled 
%to get from x to y. H returns the P(y->x)/P(x->y).
%2. When sampling directly from W(M(x)), H is the empty array []

%m: The shape(s) of which we want to find the volume are given by M(x)\leq
%m if increasing is true, and M(x)\geq m if increasing is false. m can be a
%list

%M takes in a vector x outputs the 'smallest' shape which contains x. Formally,
%M outputs \max_i g_i(x)/k_i if increasing is true, and \min_i g_i(x)/k_i if increasing is false 
%reference_volume is the volume of the largest shape 

%return_samples: if return_samples is true, we return the samples we drew along with the 
%random variables in the numerator/denominator of our sum, and which samples
%were accepted. Otherwise, we return []


%reference_volume: the volume to which we are comparing the other volumes.
%This quanity is for normalizing. 

%there are a number of options for the last two arguments:
%1. they can be unspecified, in which case the weight functions is assumed
%to be the identity
%2. Just W can be specified, in which case it is assumed to be a
%function handle and becomes the importance weights
%3. both W and nodes can be specified, in which case they are assumed to
%be numeric vectors of the same length. The weight function is
%then the trapezoidal rule for (W_in,nodes).


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator

%accepted is whether or not the proposal was accepted at each monte carlo step
function [volume,samples,numerator_rv,denominator_rv,accepted] = volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal_function,H,m,M, return_samples,reference_volume,W,nodes)
    MH=isa(H,'function_handle'); %'true' if we need a metropolis-hastings step
    %the next block of code defines the weight functions
    
    %checks for option 1 listed in the comments above
    if ~exist('W','var')%TODO consider adding checks if inputs are correct
        W_process=@M;
        %checks for options 3 in the comments above
    elseif exist('nodes','var') &&~isa(W,'function_handle')
        %if some points have infinite weight, we replace infinity values with
        %twice the maximum of the non-infinte values
        ind=W==inf;
        W(ind)=max(W(~ind))*2;
        W_process=@(y)trapezoidal_interpolation(nodes,W,M(y));
     %option 2 remains
    else
        W_process=@(y)W(M(y));
    end
    if MH
        W_importance=W_process;
    else
        W_importance=@(y)1;
        H=@(x,y,k)1;
    end
    W_estimator=@(y)1/W_process(y);
       
    [ratio,samples,numerator_rv,denominator_rv,accepted,~]  = volume_set_indicators(N,burn_in,up,increasing,x0, next,start_coordinates, proposal_function,H,m,M, return_samples,W_importance,W_estimator);   
       

    volume=ratio*reference_volume;

end