%given non uniformly spaced nodes, returns the weights for quadrature using
%the trapezoidal rule
function [weights] = non_unif_trapezoidal_weights(nodes)
n=length(nodes);
deltas=nodes(2:n)-nodes(1:(n-1));
weights=1/2*([deltas ;0]+[0; deltas]);
end

