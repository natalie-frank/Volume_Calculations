 %%inputs%%
%N is the total number of samples

%m is a numeric double

%increasing is if W is increasing or decreasing

%W is either 
%           1. a function handle
%           2. a list of postive monotonic doubles, the falue of a function of a
%function at points of m. W is assumed to be monotonic

    %%outputs%%
%Ns is a list of integers that sum to N
%Specifically, Ns puts the same weight at index i as importance sampling with N points 


function [Ns] = get_Ns_umbrella(N,m,increasing, W)
    if length(m)==1
        Ns=N;
        warning('the variable m only had one point')
    else
        
        num_points=length(m);
        if isa(W,'function_handle')
            W_list=zeros(num_points,1);
            for i=1:num_points
                W_list(i)=W(m(i));
            end
        else
            W_list=W;
        end
        if increasing
            increasing_multiplier=1;
        else
            increasing_multiplier=-1;
        end
        diff=increasing_multiplier*(W_list(2:num_points)-W_list(1:(num_points-1)))./(m(2:num_points)-m(1:(num_points-1)));
        Ns_estimate=diff/sum(diff)*N;
        Ns=round(Ns_estimate);
        Ns_difference=N-sum(Ns);
        indices=1:(num_points-1);
        if Ns_difference~=0
            if Ns_difference<0
                rounded_up=Ns>Ns_estimate;%the indices which were rounded up
                indices=indices(rounded_up);
                sign=-1;
            elseif Ns_difference>0
                rounded_down=Ns<Ns_estimate;%the indices which were rounded down
                indices=indices(rounded_down);
                sign=+1;
            end
            for i=1:abs(Ns_difference)
                index=indices(i);
                Ns(index)=Ns(index)+sign;
            end
        end
    end
end

