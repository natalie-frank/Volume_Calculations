%this is the inverse time update to f-- namely f=1/N on a log scale
%function [f_new,Hist_new,other_new] = f_update_inverse_time(f,g,N,Hist,r,dst, other,log_scale)
function [f_new,Hist_new,other_new] = f_update_inverse_time(f,g,N,Hist,r,dst, other)
    f_new=exp(1/N);
    Hist_new=Hist;
    other_new=other;
end

