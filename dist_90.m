%y is a vector of even dimension on [0,1]^d
%it is assumed that y represents disk centers on the torus [0,1]x[0,1], odd entries are
%the x-values and even entries are the y-values

%returns the minimum pairwise distance between the disk centers on
%[0,1]x[0,1]
function [d] = dist_90(y)
    if length(y)==2
        d=.5;
    else
        n=length(y)/2;
        bound=3;
        d=min(bound*eye(n)+dist_pairs_90(y),[],'all')/2;
    end
end

