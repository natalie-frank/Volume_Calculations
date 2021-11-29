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