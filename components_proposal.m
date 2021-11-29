function [y] = components_proposal(x,k,m,A)
    d=length(x);
    y=zeros(d,1);
    x_disks=x(1:(d-2));
    theta=rand(1,1)*2*pi;%TODO is this the best way to distribute theta? I suspect this is optimal because it undoes the condition number. warning: changing this will also change H
    x_move=maximum_move_disks(x_disks,k,m,theta,A); %calculating the maximum distance disk k could travel before hitting another disk, (bounded above by some constant)
    %TODO remove opt (true here)
    
    rr=x_move*(rand(1,1)); %randomly picking a distance to move
    delta=rr*[cos(theta);sin(theta)];
    y(1:(d-2))=x_disks;
    y(k)=mod(y(k)+delta,1);%moving disk k
    y(d-1)=theta;
    y(d)=x_move;    
end

