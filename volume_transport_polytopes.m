        %%input variables%%
%mu, nu: the transport polytope is given by all matrices M which satisfy  M\one=\mu, \one^T M=\nu

%N is the number of samples

%burn_in is the burn_in time

%subspace_normalization is a boolean. Let n be the dimension of mu and m
%the dimension of nu. 'subspace_normalization=true' normalizes the volume
%of the transport polytope as a subspace of R^(mn), while
%'subspace_normalization=false' calculates the R^((m-1)(n-1)) volume

%return_samples: if return_samples=true returns
%samples,numerator_rv,denominator_rv,accepted are returned

%W is the importance weight function

%node_num: 'node_num=[]' uses metropolis hastings steps in the algorithm.
%if 'node_num' is some integer, we sample from the importance distribution
%directly. However, we need to invert W numerically and we use node_num
%nodes for this.




        %%outputs%%
%volume is a vector which contains the estimate of the volume at each m

%samples is a list of our samples in the method

%numerator_rv is the realizations of the random variables in the numerator
%of our estimator

%denominator_rv is the realizations of the random variables in the
%denominator of our estimator







%TODO add 'nodes' option after W
function [volume,samples,numerator_rv,denominator_rv,accepted] = volume_transport_polytopes(mu,nu,N, burn_in,subspace_normalization,return_samples, W,node_num)
    cols=length(nu);
    [~,I_max]=max(nu);
    first_I_max=I_max(1);
    nu=[nu( (1:cols)~=first_I_max);nu(first_I_max)];
    rows=length(mu);
    [~,I_max]=max(mu);
    first_I_max=I_max(1);
    mu=[mu((1:rows)~=I_max(1)); mu(first_I_max)];
    log_vol_row_simplex=(cols-1)*sum(log(mu(1:(rows-1))))-(rows-1)*log_factorial(cols-1);
    log_vol_col_simplex=(rows-1)*sum(log(nu(1:(cols-1))))-(cols-1)*log_factorial(rows-1);
    
    %if the volume of the simplex associated with the rows is larger
    %than the simplex associated with the columns, we switch the rows and
    %columns
    if log_vol_col_simplex<log_vol_row_simplex
       temp=nu;
       dim_temp=cols;
       simplex_temp=log_vol_col_simplex;
       nu=mu;
       cols=rows;
       log_vol_col_simplex=log_vol_row_simplex;
       mu=temp;
       rows=dim_temp;
       log_vol_row_simplex=simplex_temp;
    end

    
    next=@(k)next_random(k,rows-1,cols-1); %TODO test if random or deterministic tends to do better
    start_coordinates=[1,1,1];
    samples_cell=true;
        
    if isempty(node_num)
        H=@(x,y,k)1;
        proposal_function=@(X,k)unif_simplex_proposal(X,k,mu,nu);
    else
       H=[];
       proposal_function=@(X,k)simplex_importance_proposal(X,k,mu,nu,W,node_num);
    end
    m=1;
    r=rows+1;
    c=cols+1;
    M=@(X)M_shape(X,r,c);

    increasing=true;
    reference_volume=1;
        
    X0=zeros(rows+1,cols+1);
    X0(1,1)=min(mu(1),nu(1));
    r=rows-1;
    c=cols-1;
    X0(1:r,c+1)=mu(1:r)-X0(1:r,1:c)*ones(c,1);
    X0(r+1,1:c)=ones(1,r)*X0(1:r,1:c)./nu(1:c)';
    X0(r+1,c+1)=sum(X0(1:r,1:c),'all');
    mx=max(X0(r+1,1:c));
    X0(r+2,c+2)=max(mx,(1-mu(rows)-nu(cols))/X0(r+1,c+1));
    up=true;
    
          
    [ratio,samples,numerator_rv,denominator_rv,accepted] = volume_marginal(N,burn_in,up,increasing,X0, next,start_coordinates, proposal_function,H,m,M, return_samples, reference_volume,W);

    volume=log(ratio)+log_vol_row_simplex;
    if subspace_normalization
        volume=volume+(rows-1)/2*log(cols)+(cols-1)/2*log(rows);
    end
        
        
        
    
end

function[mx]=M_shape(X,r,c)
    mx=X(r,c);
end
function [k_next]=next_random(k,rows,cols)
    %k_next=[randi(rows,1,1),randi(cols,1,1), mod(k(3),rows*cols)+1];
    k_next=[randi(rows,1,1),randi(cols,1,1), 1];
end

function [Y]=unif_simplex_proposal(X,k,mu,nu)
    Y=X;
    rows=length(mu)-1;
    cols=length(nu)-1;
    mx=X(k(1),cols+1)+X(k(1),k(2));%the maximum possible value of the coordinate we are proposing to move
    rnd=mx*rand(1,1);
    Y(k(1),k(2))=rnd;%moving the coordinate we are proposing
    if k(3)==1
        Y(1:rows,cols+1)=mu(1:rows)-Y(1:rows,1:cols)*ones(cols,1);
        Y(rows+1,1:cols)=ones(1,rows)*Y(1:rows,1:cols)./nu(1:cols)';
        Y(rows+1,cols+1)=sum(Y(1:rows,1:cols),'all');
        
    else 
        Y(k(1),cols+1)=X(k(1),cols+1)-rnd+X(k(1),k(2));%updating the row slack variable
        Y(rows+1,k(2))=X(rows+1,k(2))-rnd+X(k(1),k(2));%updating the column slack variable
        Y(rows+1,cols+1)=X(rows+1,cols+1)-X(k(1),k(2));%updating the remaining slack variable
    
    end
    %TODO there is a faster way to do this
    mx=max(Y(rows+1,1:cols));
    mx=max(mx,(1-mu(rows+1)-nu(cols+1))/Y(rows+1,cols+1));
    Y(rows+2,cols+2)=mx; 

    
end


%importance sampling from the simplex
function [Y]=simplex_importance_proposal(X,k,mu,nu,W,node_num)
    Y=X;
    rows=length(mu)-1;
    cols=length(nu)-1;
    mx=X(k(1),cols+1)+X(k(1),k(2));%the maximum possible value of the coordinate we are proposing to move
    col_sum_without_rc=X(cols+1,k(2))*nu(k(2))-X(k(1),k(2));
    all_sum_without_rc=X(rows+1,cols+1)-X(k(1),k(2));
    indices=1:cols;
    wch=indices~=k(2);
    indices=indices(wch);
    max_without_col=max(X(rows+1,indices));
    x_rc_values=linspace(0,mx,node_num);
    pdf_weights=W(M_restricted(x_rc_values,k,col_sum_without_rc,all_sum_without_rc,max_without_col,mu,nu));
    pdf_weights=make_finite(pdf_weights);

    
    
    %we now integrate by trapezoidal quadrature;
    cdf_weights=cumsum(pdf_weights);
    cdf_weights=cdf_weights-pdf_weights(1)/2;
    cdf_weights=cdf_weights-pdf_weights/2;
    cdf_weights=make_finite(cdf_weights);
    cdf_weights=cdf_weights/cdf_weights(node_num);
    [cdf_weights, ind]=unique(cdf_weights);
    x_rc_values=x_rc_values(ind);
    rnd=rand(1,1);
    new_x_rc=interp1(cdf_weights,x_rc_values,rnd,'linear');
    if ~isreal(new_x_rc)
        tt=1;
    end
    Y(k(1),k(2))=new_x_rc;%moving the coordinate we are proposing
    if k(3)==1
        Y(1:rows,cols+1)=mu(1:rows)-Y(1:rows,1:cols)*ones(cols,1);
        Y(rows+1,1:cols)=ones(1,rows)*Y(1:rows,1:cols)./nu(1:cols)';
        Y(rows+1,cols+1)=sum(Y(1:rows,1:cols),'all');
        
    else 
        Y(k(1),cols+1)=X(k(1),cols+1)-new_x_rc+X(k(1),k(2));%updating the row slack variable
        Y(rows+1,k(2))=X(rows+1,k(2))-new_x_rc+X(k(1),k(2));%updating the column slack variable
        Y(rows+1,cols+1)=X(rows+1,cols+1)-X(k(1),k(2));%updating the remaining slack variable
    
    end
    %TODO there is a faster way to do this
    mx=max(Y(rows+1,1:cols));
    mx=max(mx,(1-mu(rows+1)-nu(cols+1))/Y(rows+1,cols+1));
    Y(rows+2,cols+2)=mx; 
    
end

function[v]=make_finite(w)
    ind=w==inf | w<0;
    if all(ind)
        w=ones(node_num,1);
    elseif max(w(~ind))==0
        w(ind)=1;
    else
        w(ind)=2*max(w(~ind));
    end
    v=w;
end

function [y]=M_restricted(x_rc,k,col_sum_without_rc,all_sum_without_rc,max_without_col,mu,nu)
    m=length(nu);
    n=length(mu);
    y=zeros(length(x_rc),1);
    mult=1-nu(length(nu))-mu(length(mu));
    if m==2 &&n==2
        y=x_rc/nu(1);
    else
        if m==2
            t1=0;
            t2=0;
        else
            t1=nu(k(2))*max_without_col-col_sum_without_rc;
            t2=mult/max_without_col-all_sum_without_rc;
        end
        if n==2
            t1=0;
            t3=inf;
        else
            t3=1/2*(-(all_sum_without_rc+col_sum_without_rc)+sqrt( (all_sum_without_rc-col_sum_without_rc)^2+4*nu(k(2))*mult));
        end
        if t1<t2
            ind=x_rc<t3;
            y(ind)=mult./(all_sum_without_rc+x_rc(ind));
            y(~ind)=(col_sum_without_rc+x_rc(~ind))/nu(k(2));
        else
            ind2=x_rc<t2;
            ind1=x_rc>t1;
            ind3=~ind2 & ~ind1;
            y(ind2)=mult./(all_sum_without_rc+x_rc(ind2));
            y(ind1)=(col_sum_without_rc+x_rc(ind1))/nu(k(2));
            y(ind3)=max_without_col;
        end
        y=make_finite(y);
    end



end





%returns the log of n!
function [y]=log_factorial(n)
    v=1:n;
    y=sum(log(v));
end
