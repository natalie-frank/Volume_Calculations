
        %%inputs%%
%N is the number of samples

%burn_in is the burn_in time

%up: up is true if we compute the volume of our desired shape by comparing
%to a larger shape and false if we compar to a smaller shape

%increasing: increasing is true if our shape is given by the set of points
%which satisfy g_i(x)\leq k_i and false the shape is given by g_i(x)\geq
%k_i. Note that increasing is true if the set of shapes given by M(x)\leq
%increases with m and is false if it decreases with m

%x0 is a poin that is included in all of the shapes

%next: This is a mechanism that allows varying the coordinates in a gibbs
%sampler. Specifically, next is a function that takes the coordinates for
%which the method just make a propsal as input and outputs the next
%coordinates for which the MCMC method should make a proposal

%start_coordinates: start_coordinates is the first set of coordinates we
%use in our Gibbs sampler

%proposal_function: this function takes in three arguments (x,k,m), the current
%point and the coordinates from which to sample. for a fixed m,
%prposoal_function(x,k,m)& H(x,y,k,m) should sample uniformly from the mth
%shape


%H: H takes in 4 arguments (y,x,k,m). y is the proposal 
%move, x is the current point, and k is the coordinates that were sampled 
%to get from x to y. H returns the P(y->x)/P(x->y) for points y,x in the shape defined by m.


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


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is the samples generated by our method, returned as a cell array,
%where each entry in the cell array corresponds to a different m

%ratio_list is the list of the volume ratios which were generated by the method

%accepted is whether or not the proposal was accepted at each monte carlo
%step. Also returned as a cell array.


function [volume,samples,ratio_list,accepted] = volume_umbrella_sampling(Ns,burn_ins,up,increasing,x0, next,start_coordinates, proposal_function,H,m,M, return_samples,reference_volume)
    ratio_number=length(m)-1;
    ratio_list=zeros(ratio_number+1,1);
    if up
        add_index=0;%helps keep track of the index of ratio_list
        ratio_list(ratio_number+1)=reference_volume;
        direction='reverse';%the direction for taking the cumulative product when finding the volume
    else
        add_index=1;%helps keep track of the index of ratio_list
        ratio_list(1)=reference_volume;
        direction='forward';%the direction for taking the cumulative product when finding the volume
    end
    if ~increasing%if shapes decreasse with m, we flip all our variables so that after they will increase with m
        m=flip(m);
        burn_ins=flip(burn_ins);
        Ns=flip(Ns);
    end
    if return_samples
        samples=cell(ratio_number,1);
        accepted=cell(ratio_number,1);
    else
        samples=[];
        accepted=[];
    end
    
    initial_point=x0;
    %the reference volume is either the smallest or largest ratio
    for i=1:ratio_number
        m_sample=m(i+1);
        m_indicator=m(i);
        proposal_function_m=@(x,k)proposal_function(x,k,m_sample);
        H_m=@(x,y,k)H(x,y,k,m_sample);
        N=Ns(i);
        burn_in=burn_ins(i);
        W_importance=@(m)1;
        W_estimator=@(m)1;
        if return_samples
            [ratio_list(i+add_index),samples{i},~,~,accepted{i},initial_point] = volume_set_indicators(N,burn_in,up,increasing,initial_point, next,start_coordinates, proposal_function_m,H_m,m_indicator,M, return_samples,W_importance,W_estimator);
        else
            [ratio_list(i+add_index),~,~,~,~,initial_point] = volume_set_indicators(N,burn_in,up,increasing,initial_point, next,start_coordinates, proposal_function_m,H_m,m_indicator,M, return_samples,W_importance,W_estimator);
        end
    end
    volume=cumprod(ratio_list,direction);
    if ~increasing
        volume=flip(volume);
        samples=flip(samples);
        accepted=flip(accepted);
        ratio_list=flip(ratio_list);
    end
    
    

end