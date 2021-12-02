
        %%inputs%%
%N is the number of samples

%burn_in is the burn_in time

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

%forget_rate is a parameter which controls how quickly new information from the samples is 
%incorporated into our volume estimates

%w0 is an initial guess for the inverse volumes 


        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%z is the running volume estimate in the method

%w is the running inverse volume estimate in the method

%samples is a list of our samples in the method

%accepted is whether or not the proposal was accepted at each monte carlo step
                                
function [volume,z,w,samples,accepted] = volume_adaptive_integral_inverse_volumes(N,burn_in,increasing,x0,next,start_coordinates,proposal_function, H,m,M,return_samples,reference_volume,forget_rate, w0)
    %here we determine is we should return samples is a cell array or a
    %matrix
    %specifcially, we use a matrix if the samples are column vectors and
    %otherwise we use a cell array.
    if isnumeric(x0)
       sz=size(x0);
       if sz(2)==1
           samples_cell=false;
       end
    else
       samples_cell=true;
    end                             

    if increasing
        increasing_multiplier=1;
    else
        increasing_multiplier=-1;
    end                               
    quadrature_weights=non_uniform_trapezoidal_weights(m);
    
    
    mg=dot(w0,quadrature_weights);
    if mg~=1
        %warning(strcat('w0 numerically integrated to ',num2str(mg),' and was renormalized'));
        w0=w0/mg;
    end
    
    w=w0;%initial weight vector
    z=zeros(length(m),1);%the vector of volumes
    d=length(x0);%dimension of the problem
    if return_samples
        if samples_cell
            samples=cell(N,1);
        else
            samples=zeros(d,N);%to hold the samples we generate
        end
        accepted=zeros(N,1);%whether or not a proposal is accepted
    else
        samples=[];
        accepted=[];
    end
    current=x0;%keeping track of the current sample
    next_coordinates=start_coordinates;%keeps track of the next coordinates we want to move
    
    %dummy indicator functions
    len=length(m);
    ind_den=@(x,r) ones(len,1);
    ind_num=@(x,r) ones(len,1);


%now we take burn_in steps without including them in the calculation
    if burn_in>0
        W=@(y)quadrature_function(m,quadrature_weights,w,M(y));
        [~,samples,~,~,~,~]=volume_ratio(burn_in,0,current, next,next_coordinates, proposal_function,H,m,W,W,ind_num,ind_den,true);
        if samples_cell
            sample=samples{burn_in};
        else
            sample=samples(:,burn_in);
        end
        current=sample;
    end

    for j=1:N
        W=@(y)quadrature_function(m,quadrature_weights,w,M(y));%updating the weight function
        %W=@(y)trapezoidal_interpolation(nodes,w,dist(y));
        %we use 'find_vol_r' to generate samples, we don't care about the
        %volume outputs
        [~,~,~,~,acc,sample] = volume_ratio(1,0,current, next,next_coordinates, proposal_function,H,m,W,W,ind_num,ind_den,true);
        if return_samples
            accepted(j)=acc;
            if samples_cell
                samples{j}=sample;
            else
                samples(:,j)=sample;
            end
        end
        current=sample;
        next_coordinates=next(next_coordinates);
        
        indicator=increasing_multiplier*M(sample)<=increasing_multiplier*m;
        
        
        z_next=j/(j+1)*z+1/(j+1)*indicator/W(sample);
        %since we are about to take the inverse of z_next, we replace all zeros
        %in z_next with 1/2 of the next smallest value. This is the update of
        %our volumes z.
        z=z_next;
        index= z_next==0;
        mn=min(z_next(~index));
        z(index)=mn/2;
        
        y=1./z;%inverse of the volume
        
        %the next two lines are from the second option of the 'Discretization'
        %section of our writeup
        w_star=(1-forget_rate)*w+forget_rate*y;
        w=w_star/dot(quadrature_weights,w_star);
    
    end

%We normalize by the largest entry of z. If z is increasing with r,
%the largest entry of z should be the last one.
%Otherwise,it's the first entry.
    if increasing
        volume=z/z(length(m));
    else
        volume=z/z(1);
    end
    volume=volume*reference_volume;
end





