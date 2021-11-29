
         %%inputs%%
%N is the number of samples
%burn_in is burn in time
%x0 is a point inside the component
%m is a list of radiuses for which we want to find the volume of the component
%torus is either 60 or 90, for the 60 degree torus and the 90 degree torus
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

%if m in m_list is greater than or equal to the maximum radius at which disks x0 don't overlap, the algorithm will not run-- it will
%return 0 and leave empty the corresponding entry in samples, numerator_rv,
%denominator_rv
function [volume,samples,numerator_rv,denominator_rv,accepted] = volume_disks_components(N,burn_in,x0,m_list,torus,return_samples, W,nodes)
    p=inf;%we will be considering 2, infinity group norm
    hidden=2;
    vol_torus=1;%volume of our torus
    increasing=true;
    up=false;
    if torus==90
        M_shape=@dist_90;
        A=eye(2,2);%the matrix that maps [0,1]x[0,1] onto our torus
        fundamental=1;%volume of [0,1]x[0,1] mapped by A
    elseif torus==60
        M_shape=@dist_60;
        A=[1 1/2;%the matrix that maps [0,1]x[0,1] onto our torus
        0 sqrt(3)/2];
        fundamental=sqrt(3)/2;%volume of [0,1]x[0,1] mapped by A
    else
       error('torus must be either 60 or 90') 
    end
    proposal_function=@(x,k,m)components_proposal(x,k,m,A);
    H=@(x,y,k,m)H_components_proposal(x,y,k,A,m);
    %the next block of code processes W_in which will later be used to
    %define W
    if isa(W,'cell')
        per_m=true;
    else
        per_m=false;
    end
    

    d=length(x0);
    n=d/2;
    next= @(k) nxt(k,d);
    start_coordinates=[1;2];
    len=length(m_list);
    volume=zeros(len,1);
    
    %preallocating the output variables
    if return_samples
        denominator_rv=zeros(len,N);
        accepted=false(len,N);
        samples=cell(len,1);
    else
        denominator_rv=[];
        accepted=[];
        samples=[];
    end

      
   
    for k=1:len
        m=m_list(k);
        if per_m
            W_m=W{k};
        else
            W_m=W;
        end
        
        m_constraint=max(M_shape(x0)-m,0);
        if m_constraint<=5*eps %we use 5*eps instead of 0 to guard against numerical error
             volume(k)=0;
        else
            x0_m=x0;
            x0_m(1)=x0_m(1)+m_constraint/2;%we can't start at x0 because that point has zero weight
            reference_volume=(pi*m_constraint^2)^n*(vol_torus/fundamental)^n; 
            x0_m((d+1):(d+hidden))=zeros(hidden,1);
            proposal_function_m=@(x,k)proposal_function(x,k,m);
            H_m=@(x,y,k)H(x,y,k,m);
            M_constraint=@(x)group_norm(x(1:d),x0,p,M_shape);
            if exist('nodes','var')
                if return_samples
                    [volume(k),samp,numerator_rv,denominator_rv(k,:),accepted(k,:)]=volume_marginal(N,burn_in,up,increasing,x0_m,next,start_coordinates, proposal_function_m,H_m,m_constraint,M_constraint, return_samples,reference_volume,W_m,nodes);%finds the mean ratio ind_num/ind_den for the parameters given
                else
                    [volume(k),~,~,~,~]=volume_marginal(N,burn_in,up, increasing,x0_m,next,start_coordinates, proposal_function_m,H_m,m_constraint,M_constraint,return_samples,reference_volume, W_m,nodes); %finds the mean ratio ind_num/ind_den for the parameters given
                end
            else
                if return_samples
                    [volume(k),samp,numerator_rv,denominator_rv(k,:),accepted(k,:)]=volume_marginal(N,burn_in,up,increasing,x0_m,next,start_coordinates, proposal_function_m,H_m,m_constraint,M_constraint, return_samples,reference_volume,W_m); %finds the mean ratio ind_num/ind_den for the parameters given
                else
                    [volume(k),~,~,~,~]= volume_marginal(N,burn_in,up, increasing,x0_m,next,start_coordinates, proposal_function_m,H_m,m_constraint,M_constraint,return_samples,reference_volume, W_m);%finds the mean ratio ind_num/ind_den for the parameters given
                end
            end
            if return_samples
                samples{k}=samp(1:d,:);
            end
        end
    end
end




function [nrm]=group_norm(x1,x2,p,dist)
    n=length(x1)/2;
    y=zeros(2*n,1);
    for ii=1:n
        y(ii)=pair_dist(x1(2*ii-1: 2*ii),x2(2*ii-1: 2*ii),dist);
    end
    nrm=norm(y,p);
end

function[dst]=pair_dist(p1,p2,dist)
    dst=2*dist([p1;p2]);%should there be a factor of 2 here?
end


function [k_next]=nxt(k_begin,dim)
    ell=length(k_begin);
    k_next=k_begin+ell-1;
    k_next=mod(k_next,dim)+1;
end
