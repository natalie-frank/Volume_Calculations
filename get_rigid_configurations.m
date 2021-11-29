%n is the number of disks
%returns the rigid configurations on the 60-degree hexagon-torus with n disks and the
%corresponding radiuses
%type can be 'radius' or 'density'
%options for domain are 'tor60' or 'tor60_on_90'
%returns a matrix representing the rigid configurations, with
%a single column being a single rigid configuration
% the second return value is the raduis corresponding to each configuration
function [configurations,radiuses] = get_rigid_configurations(n,type,domain)
if ~strcmp(domain, 'hex') &&~strcmp(domain,'tor60_on_90') &&~strcmp(domain,'tor60')
    error('invalid domain')
end
if ispc
    M1=csvread('Disk_Data\Disk_Coordinates.csv',1,0);
else
    M1=csvread('Disk_Data/Disk_Coordinates.csv',1,0);
end
n1=M1(:,1);
ind1=n1==n;%all the indices that correspond to n disks
conf1=M1(:,2);%a reference number for the configuration is in the second column of the file
x=M1(:,3);%x-values are in the 3rd column of the file
y=M1(:,4);%y-values are in the 4th column of the file
conf1=conf1(ind1);%the cofiguration reference number for all critical configurations with n disks
x=x(ind1);%the x-values that correspond to n disks
y=y(ind1);%the y-values that corresond to n disks

if ispc
    M2=csvread('Disk_Data\critical_configuration_data.csv',1,0);
else
    M2=csvread('Disk_Data/critical_configuration_data.csv',1,0);
end
n2=M2(:,1);%records the number of disks in each critical configuration in the second file
ind1=n2==n;%finds the indices for which the configuration is for n disks
conf2=M2(:,2);%2nd column is configuration reference number
index=M2(:,5);%5th column is index
if strcmp(type,'radius')%radiuses recorded in 3rd column of file
    jj=3;
elseif strcmp(type, 'density')%density recorded in 4th column of file
    jj=4;
end
rads=M2(:,jj);

%extracting the relevant information just for n disks
conf2=conf2(ind1);
index=index(ind1);
rads=rads(ind1);

%finding the index zero critical points
ind2=(index==0);
conf=conf2(ind2);
rads=rads(ind2);

ind3=ismember(conf1,conf);%checks which cofigurations correspond to index 0
conf_only=conf1(ind3);
x=x(ind3);
y=y(ind3);
num=length(x)/n;
configurations=zeros(2*n,num);
radiuses=zeros(1,num);
even=2*(1:n);
odd=even-1;
if ~exist('domain','var')
    domain='hex';
end
if strcmp(domain,'tor60')||strcmp(domain,'tor60_on_90')
    v= hexagon_torus_to_60_torus([x y],'xy');
    x=v(:,1);
    y=v(:,2);
    len=length(x);
    even2=2*(1:len);
    odd2=even2-1;
    w=zeros(2*len,1);
    w(odd2)=x;
    w(even2)=y;
    v=torus_60_to_torus_90(w);
    if strcmp(domain,'tor60')
        v=torus_90_torus_tor_60(v);
    end
    x=v(odd2);
    y=v(even2);
end
for k=1:num
   configurations(odd,k)=x(((k-1)*n+1):k*n);
   configurations(even,k)=y(((k-1)*n+1):k*n);  
   config_num=conf_only((k-1)*n+1);
   ind= conf==config_num;
   radiuses(k)=rads(ind);
   
end
end



%maps the hexagon-torus to the 60 degree torus
function [v] = hexagon_torus_to_60_torus(X,format)
if strcmp(format,'xy')
    x=X(:,1);
    y=X(:,2);
elseif strcmp(format,'col')
    sz=size(X);
    n=sz(1)/2;
    even=2*(1:n);
    odd=even-1;
    x=X(odd,:);
    y=X(even,:);
    v=zeros(sz(1),sz(2));
else
    error('incorrect format string')
end
ind1= y<-1/(2*sqrt(3));
ind2= sqrt(3)*(x+1/2)-1/(2*sqrt(3))<y;


x(ind1)=x(ind1)+1/2;
y(ind1)=y(ind1)+sqrt(3)/2;
x(ind2)=x(ind2)+1;

x=x+1/2;
y=y+1/2/sqrt(3);
if strcmp(format,'xy')
    v=[x, y];
else
    v(odd,:)=x;
    v(even,:)=y;
end

end

