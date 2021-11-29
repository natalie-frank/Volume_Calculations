
%outputs P(y\to x)/P(x\to y) under the proposal
%components_proposal_original
%TODO define variables
%points in form [coordinates, theta, x_move]
function [ratio]= H_components_proposal(y,x,k,A,m)
   %if we get from x to y by moving disk k in a direction with angle theta,
    %we get from y to x by moving disk k in direction -theta
    %note P(x \to y)=1/max_dist(x,k,r,theta,A)
    %ratio=xxx_moveee/max_dist(y,k(2)/2,r,mod(pi+ttthetaaa,2*pi),A,opt);
    d=length(x);
    theta=y(d-1);
    x_move=y(d);
    y_move=maximum_move_disks(y(1:(d-2)),k,m,mod(pi+theta,2*pi),A);
    if x_move==0
        ratio=1;
    else
        ratio=x_move/y_move;
        %TODO get rid of opt (true here)
        %%%%%%%%ratio=x_move/max_dist(x,k,m,mod(pi+theta,2*pi),A,true);
    end
end

