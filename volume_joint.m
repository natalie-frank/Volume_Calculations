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


%reference_volume: the volume to which we are comparing the other volumes.
%This quanity is for normalizing. 

%W is a function handle which is the integral of the importance weights on
%the 'simulated tempering' variable. W must be increasing if the variable
%'increasing' is false, and otherwise if must be decreasing

%W_inverse is the inverse of W


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator

%accepted is whether or not the proposal was accepted at each monte carlo step
function [volume,samples,numerator_rv,denominator_rv,accepted] = volume_joint(N,burn_in,up,increasing, x0, next,start_coordinates,p, proposal_function,H,m,M, return_samples,reference_volume,W,W_inverse)

    if increasing
        increasing_multiplier=1;
    else
        increasing_multiplier=-1;
    end
      

    x0_joint=[x0;rand_shape_index(x0,M,W,W_inverse)];
    proposal_function_joint=@(x,k)joint_proposal(p,proposal_function, W,W_inverse,M, x,k);
    H_joint=@(x,y,k) H(x(1:(length(x)-1)),y(1:(length(y)-1)),k);
    M_joint=@(x)M(x(1:(length(x)-1)));

    W_estimator=@(y) 1/W(M_joint(y));
    W_importance=@(y) increasing_multiplier*M_joint(y)<= increasing_multiplier*y(length(y));
       
                                                          
    [ratio,samples,numerator_rv,denominator_rv,accepted]  = volume_set_indicators(N,burn_in,up,increasing,x0_joint, next,start_coordinates, proposal_function_joint,H_joint,m,M_joint, return_samples,W_importance,W_estimator);
    volume=ratio*reference_volume;
    
    


end



function [y]=joint_proposal(p,proposal_function,W,W_inverse,M, x,k)
    rnd=rand(1,1);
    d=length(x)-1; %dimension of input point
    y=zeros(d+1,1);
    if rnd>p
        y(1:d)=proposal_function(x(1:d),k);
        y(d+1)=x(d+1);
    else
        y(1:d)=x(1:d);
        y(d+1)=rand_shape_index(x(1:d),M,W,W_inverse);
    end

end

function [m]=rand_shape_index(x,M,W,W_inverse)
            rnd=rand(1,1);
            m=W_inverse(rnd*W(M(x)));
end