%the update to f proposed in the paper which first described the wang
%landau method
function [f_new,Hist_new,other_new] = f_update_histogram(f,g,N,Hist,r,dst,other,log_scale, tolerance)
    [lower_index,lambda,~]=interpolate_wang_landau(g,r,dst,log_scale);
    if isempty(Hist)
        Hist_new=zeros(length(r),1);
    else
        Hist_new=Hist;
    end
    Hist_new(lower_index)=Hist_new(lower_index)+lambda;
    if lower_index<length(Hist_new)
        Hist_new(lower_index+1)=Hist_new(lower_index+1)+(1-lambda);
    end
    mn=mean(Hist_new);
    if all( abs(Hist_new-mn)<(1-tolerance)*mn)
        f_new=sqrt(f);
        Hist_new=zeros(length(r),1);
    else
        f_new=f;
    end
    other_new=other;
    
    
end

