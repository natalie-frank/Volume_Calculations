%maps the 90 degree torus to the 90 degree torus
function [t] = torus_60_to_torus_90(z)
n=length(z)/2;
odd=-1+2*(1:n)';
even=2*(1:n);
x=z(odd);
y=z(even);
x_new=x-1/sqrt(3)*y;
y_new=2/sqrt(3)*y;
t=zeros(2*n,1);
t(odd)=x_new;
t(even)=y_new;
t=mod(t,1);
end

