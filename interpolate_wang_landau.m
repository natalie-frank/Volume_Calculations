%since we are applying wang-landau to a continuous system, a data point
%typically falls bewteen two levels of m in our discretization
%this function figures out what weight to assign to those two levels of m 

function [index, lambda, value]=interpolate_wang_landau(g,r,x,log_scale)
    m=length(g);
    index=find(r<x,1,'last'); %we look for an index that satisfies r(index)<=x<r(index+1)
    if isempty(index)
       if x~=r(1)
            warning("proposed a radius below the minimum of inputed radiuses")
       end
        value=g(1);%since we can't evaluate g at x, we instead use the value at the next closest point
        lambda=1;
        index=1;
    elseif index==m
        if x~=r(m)
            warning("proposed a radius above the maximum of inputed radiuses")
        end
        value=g(m);%since we can't evaluate g at x, we instead use the value at the next closest point
        lambda=0;
    else
        lambda=(r(index+1)-x)/(r(index+1)-r(index));%for this lambda, r_proposal=lambda *r(lower_index)+(1-lambda)*r(lower_index+1)
        if log_scale
            %we split into two cases to make sure that the operations in
            %the calculation are well-conditioned
            if log(lambda)+g(index)> log(1-lambda)+g(index+1)
                value=log(lambda)+g(index)+log(1+(1-lambda)/lambda*exp(g(index+1)-g(index)));
            else
                value=log(1-lambda)+g(index+1)+log(1+lambda/(1-lambda)*exp(g(index)-g(index+1)));
            end
        else
            value=lambda*g(index)+(1-lambda)*g(index+1);%we interpolate the value of g at x
        end
        
    end

end
