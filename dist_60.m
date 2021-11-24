%y is a vector of even dimension on [0,1]^d
%it is assumed that y represents disk centers on the torus [0,1]x[0,1], odd entries are
%the x-values and even entries are the y-values

%returns the minimum pairwise distance between these disk centers assuming
%that the 90 degree torus [0,1] x[0,1] is mapped to the 60 degree torus via
%the map A=[1 1/2 
%           0 sqrt 3 /2]
function [d] = dist_60(y)
    if length(y)==2
        d=sqrt(3)/4;
    else
        bound=3;
        n=length(y)/2;
        d=min(bound*eye(n)+dist_pairs_60_on_90(y),[],'all')/2;
    end
end

