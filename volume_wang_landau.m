
                %%inputs%%
%N is the number of samples

%increasing is a boolean. If it's true, then we assume that the volume
%increases with m and decreases otherwise

%x0 is is a starting point for the monte carlo method

%start_coordinates is the first coordinate(s) we move with the proposal

%next: This is a mechanism that allows varying the coordinates in a gibbs
%sampler. Specifically, next is a function that takes the coordinates for
%which the method just make a propsal as input and outputs the next
%coordinates for which the MCMC method should make a proposal

%proposal_function: this function takes in two arguments (x,k), the current
%point and the coordinates from which to sample. it outputs the proposal
%move in the MCMC method

%H: H takes in 3 arguments (y,x,k). y is the proposal move, x is the
%current point, and k is the coordinates that were sampled to get from k to
%y. H returns the P(y->x)/P(x->y).

%f_update is a function that takes in (f,g,N,Hist,r,dst,other) It outputs
%an update to f, Hist, other (TODO improve description here)

%m discretizes the energy space

%largest_volume is for normalizing. We assume we know the volume of the largest shape 

%samples_return is true if we want to return the samples and if samples were accepted, false otherwise

%log_scale_calculations indicates if we should do our calculations on a log
%scale.

        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m
%samples is a list of our samples in the method
%accepted is whether or not the proposal was accepted at each monte carlo step
%g is the estimated derivative of the volumes
function [volume,samples,accepted,g] = volume_wang_landau(N,increasing, x0,start_coordinates,next,proposal_function,H,f_update,m,M,largest_volume, samples_return,log_scale_calculations)
    d=length(x0);%dimension of the problem
    len=length(m);%number of points in the discretization
    if samples_return
        samples=zeros(d,N);%will hold samples
        accepted=false(N,1);%will hold if samples were accepted or rejected
    else
        samples=[];
        accepted=[];
    end
    
    %initial variables for input to f_update
    Hist=[];
    other=[];
    
    if log_scale_calculations 
        g=zeros(len,1);%log of 1 is zero
    else
        g=ones(len,1);%will hold the derivative of the volume
    end
    f=exp(1);%we start the multipication factor of wang-landau at e
    k=start_coordinates;%first coordinates to modify in the proposal
    current=x0;%starting point for the markov chain
    m_current=M(current);%radius of current point

       
    for i=1:N
        %g_curr is the value of g at derivative the current point (recall we start g as the all ones vector)
        %we next find the index for which m(i)\leq M(current)< m(i+1), and the
        %lambda for which m_current=lambda*m(i)+(1-lambda)*m(i+1)
        [current_lower_index,lambda_current,g_current]=interpolate_wang_landau(g,m,M(current),log_scale_calculations);%after updating g, we need to recompute g_curr
        proposal=proposal_function(current,k);%generating the proposal
        m_proposal=M(proposal);%the distance of the proposal
        [proposal_lower_index,lambda_proposal,g_proposal]=interpolate_wang_landau(g,m,m_proposal,log_scale_calculations);%find index for which r(i)<=r_proposal<r(i+1), lambda for which r_proposal=lambda*r(i)+(1-lambda)*r(i+1), the value of g at r_proposal
        %accept reject step
        if log_scale_calculations
            ratio=exp(g_current-g_proposal);
        else
            ratio=g_current/g_proposal;
        end
        p_accept=min(ratio*H(proposal,current,k),1);
        rnd=rand(1,1);
        if rnd<=p_accept
           current=proposal;
           current_lower_index=proposal_lower_index;
           lambda_current=lambda_proposal;
           m_current=m_proposal;
           if samples_return
               accepted(i)=true;%record that proposal was accepted
           end
        end
        if samples_return
            samples(:,i)=current;%record the sample
        end
        
       %update the value of g
       if log_scale_calculations
           g(current_lower_index)=g(current_lower_index)+log(f)*lambda_current;
       else
           g(current_lower_index)=g(current_lower_index)*f^lambda_current;
       end
       if(current_lower_index<len)
           if log_scale_calculations
               g(current_lower_index+1)=g(current_lower_index+1)+log(f)*(1-lambda_current);
           else
               g(current_lower_index+1)=g(current_lower_index+1)*f^(1-lambda_current);
           end
       end

       [f,Hist,other]=f_update(f,g,i,Hist,m,m_current,other); %update the value of f
       k=next(k);%updating the next coordinates for the proposal
    end    
    weights=non_uniform_trapezoidal_weights(m);%getting the quadrature weights
    if log_scale_calculations
        log_weighted_g=g+log(weights);%'multiplying' by quadrature weights
        
        %since se re-normalize later anyways, we can at a constant before
        %applying exp() so that the values of exp(weighted_g) lie in a
        %reasonable range. we pick a constant c for which -(min log_weighted_g
        %+c)=max(log_weighted_g+c)
        %c=-(min(log_weighted_g)+max(log_weighted_g))/2;
        
        c=-max(log_weighted_g);
        log_weighted_g=log_weighted_g+c;
        weighted_g=exp(log_weighted_g);
    else
        weighted_g=g.*weights;
        c=1/max(weighted_g);
        weighted_g=c*weighted_g;
    end
    
    if increasing
        volume=cumsum(weighted_g);
    else
        volume=flip(cumsum(flip(weighted_g)));
    end
    
    %normalizing
    if increasing
        volume=volume/volume(len)*largest_volume;
    else
        volume=volume/volume(1)*largest_volume;
    end
    
    if log_scale_calculations
        g=g-max(g);
        g=exp(g);
    end
end


