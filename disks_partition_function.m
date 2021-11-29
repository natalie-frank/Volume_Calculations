         %%inputs%%
%N is the number of samples
%n is the number of disks
%m is the list of radiuses at which we want to compute the volume
%m is a list of radiuses for which we want to find the volume of the component
%torus is either 60 or 90, for the 60 degree torus and the 90 degree torus
%step_size is the maximum distance we move a disk at every MCMC step
%return_samples is true if the user wants to returns samples,
%denominator_rv, or accepted

%there are three options for the last two arguments:
%1. Just W can be specified, in which case it is assumed to be a
%function handle and becomes the integral of the weight function
%2. both W and nodes can be numeric vectors, in which case they are assumed to
%be numeric vectors of the same length. The weight function is
%then the trapezoidal rule for (W_in,nodes).
%3. nodes can be a numeric vector kx1 vector, and then W can be a cell array 
%of length length(m) of kx1 numeric vectors. Then for each element of m,
%say m(i), the weight function is trapezoidal interpolation of the vectors (W_in{i}, nodes) 



        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator

%accepted is whether or not the proposal was accepted at each monte carlo step

%if m in m_list is the radius of x0, the algorithm will not run-- it will
%return 0 and leave empty the corresponding entry in samples, numerator_rv,
%denominator_rv
function [volume,samples,numerator_rv,denominator_rv,accepted] = disks_partition_function(N,n,m,torus,step_size,return_samples,W)
    if isnumeric(W)
       if length(W)~=length(m)
           error('if W is a numeric vector, it must have the same length as m');
       end
    end
    d=2*n;
    burn_in=0;
    up=true;
    increasing=false;
    x0=rand(d,1);
    next=@(k)nxt(k,d);
    start_coordinates=[1;2];
    proposal=@(x,k)step_proposal(x,k,step_size);
    H=@(x,y,k)1;
    if torus==60
        M=@dist_60;
        reference_volume=(sqrt(3)/2)^n;
    elseif torus==90
        M=@dist_90;
        reference_volume=1;
    else
        error('must specify either the 90 degree or 60 degree torus')
    end
    if return_samples 
        if isnumeric(W)
            [volume,samples,numerator_rv,denominator_rv,accepted]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,H,m,M, return_samples,reference_volume,W,m);
        else
            [volume,samples,numerator_rv,denominator_rv,accepted]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,H,m,M, return_samples,reference_volume,W);
        end
    else
        if isnumeric(W)
            [volume,~,~,~,~]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,H,m,M, return_samples,reference_volume,W,m);
        else
            [volume,~,~,~,~]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,H,m,M, return_samples,reference_volume,W);
        end
        samples=[];
        numerator_rv=[];
        denominator_rv=[];
        accepted=[];
    end
end

function [next_coordinates]=nxt(coordinates,dim)
    ell=length(coordinates);
    next_coordinates=coordinates+ell-1;
    next_coordinates=mod(next_coordinates,dim)+1;
end

function [y]=step_proposal(x,k,step_size)
    ell=length(k);
    delta=step_size*(2*rand(ell,1)-1);
    y=x;
    y(k)=y(k)+delta;
    y=mod(y,1);
end