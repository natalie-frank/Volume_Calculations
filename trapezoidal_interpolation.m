%Uses the trapezoidal rule to interpolate function values
%nodes_old is the points where we have the function values. Listed in
%increasing order
%f_old is the function values at these points
%nodes_new is points where we want the new function values
%f_new is the interpolated function values as a column vector

%TODO: add warnings if below minimum or above maximum of nodes_old
function [f_new] = trapezoidal_interpolation(nodes_old, f_old, nodes_new)
m=length(nodes_old);
n=length(nodes_new);
f_new=zeros(n,1);
old=1;%current index for old variable
new=1;%current index for new variable
 
    

while new<=length(nodes_new)&&nodes_new(new)<=nodes_old(old)
    %for all x-values we want to interpolate that are below the given range
    %we set the function equal to the left endpoint
    f_new(new)=f_old(old);
    new=new+1;
end


while old<=m-1 && new<=n
    while old<m-1 && (nodes_new(new)>nodes_old(old+1)||nodes_old(old)==nodes_old(old+1))%this finds the two consecutive points but nonequal in nodes_old that lie on either side of new
        old=old+1;
    end
    t=(nodes_new(new)-nodes_old(old))/(nodes_old(old+1)-nodes_old(old));%perform the interpolation
    f_new(new)=(1-t)*f_old(old)+t*f_old(old+1);
    new=new+1;%increment new
end
if new<n%for all x-values we want to interpolate above the given range we set to the right endpoint
    for i=(new+1):n
        f_new(i)=f_old(m);
    end
end

end

%below are outlined two more succinct way to do this for a single point in
%nodes_new
% function [val]=linear_interpolation(nodes,W_in,y)
%         if length(nodes)~=length(W_in)
%             error("'nodes' and 'W_in' must have the same length")
%         end
%         if y>=max(nodes)
%             val=W_in(length(nodes));
%         else
%             m=find(max(nodes<=y));
%             delta=nodes(m+1)-nodes(m);
%             p=(nodes(m+1)-y)/delta;
%             val=p*W_in(m)+(1-p)*W_in(m+1);
%         end
% 
% 
% 
% end



%         lower_index=find(r<r_proposal,1,'last'); %we look for an index that satisfies r(index)<=r<r(index+1)
%         if isempty(lower_index)
%             warning("proposed a radius below the minimum of inputed radiuses")
%             g_proposal=g(1);%since we can't evaluate g at r_proposal, we instead use the value at the next closest point 
%         elseif lower_index==m 
%             if r_proposal~=r(m)
%                 warning("proposed a radius above the maximum of inputed radiuses")
%             end
%             g_proposal=g(m);%since we can't evaluate g at r_proposal, we instead use the value at the next closest point 
%         else
%             lambda=(r(lower_index+1)-r_prop)/(r(lower_index+1)-r(lower_index));%for this lambda, r_proposal=lambda *r(lower_index)+(1-lambda)*r(lower_index+1)
%             g_proposal=lambda*g(lower_index)+(1-lambda)*g(lower_index+1);%we interpolate the value of g at r_proposal
%         end
