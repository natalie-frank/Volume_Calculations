
        %%inputs%%
%N is the number of samples

%burn_in is the burn_in time

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


%W_importance is a function that takes in a sample x and outputs the
%importance weight at x

%indicator_numerator is a function that takes (x,m) as inputs and outputs a
%boolean. 
%The other methods in this project assume that it outputs true if either 
%1) the shape given by m is the larger shape in the comparison of the two 
%volumes (i.e. 'up' in volume_marginal/volume_joint is false)  
%2) x is in the shape indexed by m

%indicator_denominator is a function that takes (x,m) as inputs and outputs a
%boolean. 
%The other methods in this project assume that it outputs true if either 
%1) the shape given by m is the smaller shape in the comparison of the two 
%volumes (i.e. 'up' in volume_marginal/volume_joint is true) 
%2) x is in the shape indexed by m

%W_estimator weights our samples when calculating the ratio of our two volumes.
%Explicitly, in the 'volume_ratio' function, the volume ratio is esimated 
%as \sum_{i=1}^N indicator_numerator(x,m)* W_estimator(x) / \sum_{i=1}^N
%indicator_denominator(x,m) W_estimator(x)
%('indicator_numerator' and 'indicator_denominator' are two other inputs to
%volume_ratio)


%return_samples: if return_samples is true, we return the samples we drew along with the 
%random variables in the numerator/denominator of our sum, and which samples
%were accepted. Otherwise, we return []


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator

%accepted is whether or not the proposal was accepted at each monte carlo step
function [volume,samples,numerator_rv,denominator_rv,accepted,last_sample] = volume_ratio(N,burn_in,x0, next,start_coordinates,proposal_function,H,m,W_importance,W_estimator, numerator_indicator,denominator_indicator,return_samples)
    %here we determine is we should return samples is a cell array or a
    %matrix
    %specifcially, we use a matrix if the samples are column vectors and
    %otherwise we use a cell array.
    if isnumeric(x0)
       sz=size(x0);
       if sz(2)==1
           samples_cell=false;
       else
           samples_cell=true;
       end
    else
       samples_cell=true;
    end
    MM=length(m);%number of points in discretization of m
    %if we want the samples and acceptance probability
    if return_samples
        if samples_cell
            samples=cell(N,1);
        else
            d=length(x0);%dimension of the input space
            samples=zeros(d,N);
        end
        accepted=false(N,1);
        numerator_rv=zeros(MM,N);%to hold the indicator function for the numerator evaluated at our samples
        denominator_rv=zeros(MM,N);%to hold the indicator function for the denominator evaluated at our samples
    else
        samples=[];
        numerator_rv=[];
        denominator_rv=[];
        accepted=[];
        numerator_sum=zeros(MM,1);
        denominator_sum=zeros(MM,1);
    end
    
    next_coordinates=start_coordinates;%choosing the starting coordinates to move
    current=x0;%the point we are currently at 
    W_current=W_importance(x0);% the weight of the current point
 
    for index=1:(N+burn_in)
        proposal=proposal_function(current,next_coordinates);%generating the proposal
        W_proposal=W_importance(proposal);%Weight of the proposal
        HH=H(proposal,current,next_coordinates);
        acceptance_probability=min(1,W_proposal/W_current*HH);%acceptance probability
        rnd=rand(1,1);
        
        if rnd<acceptance_probability%if the proposal is accepted
            current=proposal;
            W_current=W_proposal;
            if return_samples &&index>burn_in
                accepted(index-burn_in)=true;
            end
        end
        if return_samples && index>burn_in%recording the sample
            if samples_cell
                samples{index-burn_in}=current;
            else
                samples(:,index-burn_in)=current;
            end
        end
        if W_current==0
            error('reached zero weight point')
        end
        %saving the values of the numerator and denominator indicator
        if index>burn_in
           if return_samples
               numerator_rv(:,index-burn_in)=numerator_indicator(current,m)*W_estimator(current);
               denominator_rv(:,index-burn_in)=denominator_indicator(current,m)*W_estimator(current);
           else
               numerator_sum=numerator_sum+ numerator_indicator(current,m)*W_estimator(current);
               denominator_sum=denominator_sum+ denominator_indicator(current,m)*W_estimator(current);
           end
       end    
       next_coordinates=next(next_coordinates);%choosing the next disk to move
    end
    %finding the ratio
    if return_samples
        volume=sum(numerator_rv,2)./sum(denominator_rv,2);
    else
        volume=numerator_sum./denominator_sum;
    end
    last_sample=current;
end