%one paper recommended starting with inverse time updates and then changing
%to the proposal suggested in the wang landau paper
function [f_new,Hist_new,other_new] = f_update_optimal(f,g,N,Hist,r,dst,other,log_scale, tolerance)
    if isempty(other)
        [f_new,Hist_new,other_new]=f_update_histogram(f,g,N,Hist,r,dst, other,log_scale,tolerance);
        if f_new<exp(1/N)
            other_new=true;
        end
    else
        [f_new,Hist_new,other_new] = f_update_inverse_time(f,g,N,Hist,r,dst,other);
    end
end
