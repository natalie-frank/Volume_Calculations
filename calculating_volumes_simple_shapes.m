rng('shuffle');
%M-functions for methods
level=.3;%point of discontinuity in the derivative in the second two shapes below
sq_handle=@(x)square(x);%this outputs |x-c|_infty where c=.5[1,1]
rect_square_handle=@(x)rect_then_square_cont(x,level);% let c=.5[1.1]. If |x(1)-c(1)|<=level, it outputs |x(2)-c(2)| and otherwise this outputs |x-c|_\infty
rect_square_nonlinear_handle=@(x)rect_then_square_cont_nonlinear(x,level);%similar to the previous function, but this is cooked up to be nonlinear for |x(1)-c(1)| between 0 and level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%the list of M-functions for which we want plots


%M_function_list={sq_handle,rect_square_handle, rect_square_nonlinear_handle, @dist_90, @dist_60};
%shape_names={'square', 'rectangle then square', 'rectangle then square nonlinear','disks90','disks60'};%gives each shape a name. Very important to make sure order of names is same order as order of function handles!
%M_function_list={sq_handle};
%shape_names={'square'};
%M_function_list={@dist_90};
%shape_names={'disks90'};
M_function_list={@dist_60};
shape_names={'disks60'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters for plots
volume_vs_r=true; %do you want to see the plot of volume vs. r?
plot_sample_coordinates=true; %do you want to see the plot of the first two coordinates of the samples?
plot_histograms=true; %do you want to see a histogram of the rs?
rejection_rates=true; %do you want to see a histogram of rejection rates?
fig_numbering=2; %matlab starts numbering figures at this number
save_figures=true;%if we want to save the figures
folder='plots_and_data\Simple_Volumes_Figures';


%the following variables determine which methods we will be comparing
uniform_method=true; %uniform W
inverse_derivative_method=true; %weight function proportional to inverse derivative
antiderivative_inverse_volume_method=true; %weight function proportional to the antiderivative of the inverse weights
adaptive_antiderivative_inverse_volume_method=true;%adaptive method converging to weights proportional to the antiderivative of the inverse weights
wang_landau_method=true;%if we want to include the wang-landau method
joint_method=true;%if we want to include the method that samples in the joint distribution


%parameters for methods-- one we may want to change
d=4;%the dimension
num_points=51;%number of points in our discretization
N=1000;%number of samples
burn_in=100;%burn_in_time
reps=10;%number of repetitions
next=@(k)nxt(k,d);%picking the next coordinates to move in gibbs samplers
start_coordinates=[1,2];%the first coordinates to move in gibbs samplers
step_size=1; %note--optimal for wang-landau in 2d seems to be .1
forget_rate_anti_adaptive=10^-8*ones(length(M_function_list),1);%rate parameter for the inverse volumes method
p=0.1;%rate of changing r in joint distribution.
fixed_bins=true;%if we want a fixed number of bins in our histograms
num_bins=15;% the number of bins we want in our histograms


%most parameters you want to change should be above this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters for methods that we probably don't want to change around
proposal=@(x,k) step_proposal(x,k,step_size);%the proposal function

%some extra parameters for the Wang-Landau method
WL_dist_discrete=false;%if we want to use the discrete version of the WL method
log_scale_calculations=true;%log scale calculations in the wang-landau method
tolerance=.2; %the some of the wang landau updates use histograms of M-values. This tolerance controls how close the histogram should be to uniform before updating f 
f_update=@(f,g,N,Hist,r,dst,other)f_update_optimal(f,g,N,Hist,r,dst,other,log_scale_calculations,tolerance);
%f_update=@(f,g,N,Hist,r,dst,log_scale,other)f_update_inverse_time(f,g,N,Hist,r,dst,log_scale,other);


%weight parameters for our methods
W_uniform=ones(num_points,1);%uniform weights
w0=ones(num_points,1); %starting point for adaptive weights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%labeling our methods
tags={};
title_tags={};
if uniform_method
   len=length(tags);
   tags{len+1}='uniform';
   title_tags{len+1}='uniform weights';
end
if inverse_derivative_method
   len=length(tags);
   tags{len+1}='inverse_derivative';
   title_tags{len+1}='inverse derivative weights';
end
if antiderivative_inverse_volume_method
   len=length(tags);
   tags{len+1}='anti_inverse_volume';
   title_tags{len+1}='integral inverse volume weights';
end
if adaptive_antiderivative_inverse_volume_method
   len=length(tags);
   tags{len+1}='adaptive_anti_inverse_volume';
   title_tags{len+1}='adaptive integral inverse volume';
end
if wang_landau_method
    len=length(tags);
    tags{len+1}='wang_landau';
    title_tags{len+1}='Wang-Landau method';
end
if joint_method
    len=length(tags);
    tags{len+1}='joint';
    title_tags{len+1}='joint distribution';
end
    
colors={'k','b','m','r','c','y'}; %plot each method in a diffferent color





%initializing matricies that will store volumes
%for k=1:length(tags)
   if uniform_method
        volume_uniform_list=zeros(num_points,reps);
   end
   if antiderivative_inverse_volume_method
       volume_anti_inverse_volume_list=zeros(num_points,reps);
   end
   if inverse_derivative_method
       volume_inverse_derivative_list=zeros(num_points,reps);
   end
   if adaptive_antiderivative_inverse_volume_method
       volume_adaptive_anti_inverse_volume_list=zeros(num_points,reps);
   end
   if wang_landau_method
       volume_wang_landau_list=zeros(num_points,reps);
   end
   if joint_method
       volume_joint_list=zeros(num_points,reps);
   end
%end

%initializing matricies to store samples
if plot_sample_coordinates || plot_histograms || rejection_rates
   if uniform_method
        accepted_uniform=zeros(d,reps*N);
   end
   if antiderivative_inverse_volume_method
        samples_anti_inverse_volume=zeros(d,reps*N);
   end
   if inverse_derivative_method
       samples_inverse_derivative=zeros(d,reps*N);
   end
   if adaptive_antiderivative_inverse_volume_method
       samples_adaptive_anti_inverse_volume=zeros(d,reps*N);
   end
   if wang_landau_method
       volume_wang_landau_list=zeros(num_points,reps);
   end
   if joint_method
       volume_joint_list=zeros(num_points,reps);
   end
    if wang_landau_method
        samples_wang_landau=zeros(d,reps*N);
    end
    if joint_method
        samples_joint=zeros(d+1,reps*N);
    end
end

%intializing vectors that will store booleans representing if a proposal
%was accepted or not
if rejection_rates
   if uniform_method
        accepted_uniform=false(N*reps,1);
   end
   if antiderivative_inverse_volume_method
        accepted_anti_inverse_volume=false(N*reps,1);
   end
   if inverse_derivative_method
      accepted_inverse_derivative=false(N*reps,1);
   end
   if adaptive_antiderivative_inverse_volume_method
      accepted_adaptive_anti_inverse_volume=false(N*reps,1);
   end
   if wang_landau_method
       volume_wang_landau_list=zeros(num_points,reps);
   end 
   if wang_landau_method
       accepted_wang_landau=false(N*reps,1);
   end
   if joint_method
       accepted_joint=false(N*reps,1);
   end
end
if ~exist(folder, 'dir')
    mkdir(folder)
end

for j=1:length(M_function_list)
   dist=M_function_list{j};
   exact_volumes=zeros(num_points,1);
   volume_derivative=zeros(num_points,1);
   %variables unique for each shape
   if isequal(dist,sq_handle)
       reference_volume=1;%reference volume
      r=linspace(0,0.5,num_points)'; %discretizing radiuses
      up=true;
      increasing=true;
      exact_volumes=(2*r).^d;
      inverse_volume_weights=r.^(-d+1);
      volume_derivative=2*d*(2*r).^(d-1);
      R_lim=[0,1/2];
   elseif isequal(dist,rect_square_handle)
       reference_volume=1;%reference volume
       r=linspace(0,0.5,num_points)'; %discretizing radiuses
       up=true;
       increasing=true;
        ind=r<=(level);
        exact_volumes(ind)=(2*level).*(2*(r(ind))).^(d-1);
        volume_derivative(ind)=(2*level)*2*(d-1)*(2*r(ind)).^(d-2);
        exact_volumes(~ind)=(2*(r(~ind))).^d;
        volume_derivative(~ind)=2*d*(2*r(~ind)).^(d-1);
        R_lim=[0,1/2];
   elseif isequal(dist,rect_square_nonlinear_handle)
       reference_volume=1;%reference volume
       r=linspace(0,0.5,num_points)'; %discretizing radiuses
       up=true;
       increasing=true;
       k=1-level;
       r_2=r.^2+k*r;
       vol_rect=(2*r).^(d-1).*(2*r_2);
       volume_rectangle_derivative=2^d*r.^(d-2).*((d-1)*r_2+r.*(2*r+k));
       vol_sq=(2*r).^d;
       volume_square_derivative=2^d*d*r.^(d-1);
       ind=r<=level;
       exact_volumes(ind)=vol_rect(ind);
       exact_volumes(~ind)=vol_sq(~ind);
       exact_volumes=min(exact_volumes,1);
       volume_derivative(ind)=volume_rectangle_derivative(ind);
       volume_derivative(~ind)=volume_square_derivative(~ind);
       R_lim=[0,1/2];
   elseif isequal(dist, @dist_90)
       reference_volume=1;%reference volume
       if d~=4
           error('for disks in a torus, exact volumes only known for 2 disks')
       end
       up=true;
       increasing=false;
       r_max=sqrt(2)/4;
       r=linspace(0,r_max,num_points)';
       ind=r>1/4;
       exact_volumes(~ind)=1-4*pi*r(~ind).^2;
       exact_volumes(ind)=1-4*pi*r(ind).^2+16*r(ind).^2 .*acos(1./(4*r(ind)))-sqrt(16*r(ind).^2-1);
       volume_derivative= 8*pi*r;
       r_above=r(ind);
       volume_derivative(ind)= volume_derivative(ind)-32*r_above.*acos(1./(4*r_above))-4./sqrt(1-1./(4*r_above).^2)+16*r_above./sqrt(16*r_above.^2-1);
       volume_derivative=max(volume_derivative,0);
       R_lim=[0,r_max];
   elseif isequal(dist,@dist_60)
       reference_volume=sqrt(3)/2;%reference volume
        if d~=4
           error('for disks in a torus, exact volumes only known for 2 disks')
       end
       up=true;
       increasing=false;
       [xs,rads] = get_rigid_configurations(d/2,'radius','tor60_on_90');
       x0=xs(:,1);
       r_max=max(rads);
       r=linspace(0,r_max,num_points)';
       ind=r>1/4;
       exact_volumes(~ind)=sqrt(3)/2-4*pi*r(~ind).^2;
       exact_volumes(ind)=sqrt(3)/2-4*pi*r(ind).^2+3/2*(16*r(ind).^2 .*acos(1./(4*r(ind)))-sqrt(16*r(ind).^2-1));
       volume_derivative= 8*pi*r;
       r_above=r(ind);
       volume_derivative(ind)=volume_derivative(ind)+3/2*(-32*r_above.*acos(1./(4*r_above))-4./sqrt(1-1./(4*r_above).^2)+16*r_above./sqrt(16*r_above.^2-1));
       R_lim=[0,r_max];
   end
   if WL_dist_discrete
        wdist=@(x)discretize(dist(x),r);
   else
        wdist=dist;
   end
   inverse_derivative=1./volume_derivative;
   q_weights=non_uniform_trapezoidal_weights(r);
   inv_vol=1./exact_volumes;
   indd=inv_vol==inf;
   inv_vol(indd)=max(inv_vol(~indd))*2;
   if increasing
      inverse_volume_weights=flip(cumsum(q_weights./flip(inv_vol)));
   else
      inverse_volume_weights=cumsum(q_weights./inv_vol); 
   end
   if joint_method
       if increasing
           m_lim=R_lim(2);
       else
           m_lim=R_lim(1);
       end
       W_joint= @(r) 1; 
       W_joint_inverse=@(r) m_lim;
   end
   
   %looping over methods
   for i=1:reps
        x0=rand(d,1);
        if ~plot_sample_coordinates && ~plot_histograms && ~rejection_rates
            if uniform_method
                volume_uniform_list(:,i)=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,W_uniform,r);
            end
            if antiderivative_inverse_volume_method
                volume_anti_inverse_volume_list(:,i)=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,inverse_volume_weights,r);
            end
            if inverse_derivative_method
                volume_inverse_derivative_list(:,i)=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,inverse_derivative,r);
            end
            if adaptive_antiderivative_inverse_volume_method
                volume_adaptive_anti_inverse_volume_list(:,i)=volume_adaptive_integral_inverse_volumes(N,burn_in,increasing,x0,next, start_coordinates,  proposal, @H,r,dist,plot_sample_coordinates,reference_volume,forget_rate_anti_adaptive(j), w0);
            end
            if wang_landau_method
                volume_wang_landau_list(:,i)=volume_wang_landau(N,increasing,x0,start_coordinates,next,proposal,@H,f_update,r,wdist,reference_volume, plot_sample_coordinates,log_scale_calculations);
            end
            if joint_method
                 volume_joint_list(:,i)=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,W_joint,W_joint_inverse);  
            end
        else
            if uniform_method
                [volume_uniform_list(:,i),samples_uniform(:,((i-1)*N+1):(i*N)),~,~,accepted_uniform(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,W_uniform,r);
            end
            if antiderivative_inverse_volume_method
                [volume_anti_inverse_volume_list(:,i),samples_anti_inverse_volume(:,((i-1)*N+1):(i*N)),~,~,accepted_anti_inverse_volume(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,inverse_volume_weights,r);
            end
            if inverse_derivative_method
                [volume_inverse_derivative_list(:,i),samples_inverse_derivative(:,((i-1)*N+1):(i*N)),~,~,accepted_inverse_derivative(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,inverse_derivative,r);
            end
            if adaptive_antiderivative_inverse_volume_method
                [volume_adaptive_anti_inverse_volume_list(:,i),~,~,samples_adaptive_anti_inverse_volume(:,((i-1)*N+1):(i*N)),accepted_adaptive_anti_inverse_volume(((i-1)*N+1):(i*N))]=volume_adaptive_integral_inverse_volumes(N,burn_in,increasing, x0,next, start_coordinates,  proposal, @H,r,dist,true,reference_volume,forget_rate_anti_adaptive(j), w0);                                                                                                         
            end
            if wang_landau_method
                [volume_wang_landau_list(:,i),samples_wang_landau(:,((i-1)*N+1):(i*N)), accepted_wang_landau(((i-1)*N+1):(i*N)),~]=volume_wang_landau(N,increasing,x0,start_coordinates,next,proposal,@H,f_update,r,wdist,reference_volume,true,log_scale_calculations);
            end
            if joint_method
                [volume_joint_list(:,i),samples_joint(:,((i-1)*N+1):(i*N)),~,~, accepted_joint(((i-1)*N+1):(i*N))]=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist,true,reference_volume,W_joint,W_joint_inverse);  
            end
        end
   end
    
   %averaging
   means=zeros(num_points, length(tags));
   sds=zeros(num_points,length(tags));
   for k=1:length(tags)
      mat=eval(strcat('volume_',tags{k},'_list'));
      means(:,k)=mean(mat,2);
      sds(:,k)=sqrt(var(mat,1,2)+(means(:,k)-exact_volumes).^2);
      %sds(:,k)=sqrt(var(mat,1,2));
   end
   
   %ploting volume vs. r
   if volume_vs_r
      fig=figure(fig_numbering);
      plot_handles=gobjects(length(tags)+1,1);
      plot_handles(1)=plot(r, exact_volumes,'Color','green');
      hold on;
      for k=1:length(tags)
         plot_handles(k+1)=plot(r,means(:,k),colors{k});
         hold on;
         plot(r,means(:,k)-2*sds(:,k),strcat(colors{k},'--'),r,means(:,k)+2*sds(:,k),strcat(colors{k},'--'));
         hold on;
      end
        if increasing
            loc='Southeast';
        else
            loc='Southwest';
        end
        legend_names=cell(length(tags)+1,1);
        legend_names{1}='true volume';
        legend_names(2:length(tags)+1)=title_tags;
        legend(plot_handles,legend_names,'Location',loc);
        xlabel('r');
        ylabel('volume');
        title(strcat('volume vs. r--',shape_names{j}));
        hold off;
          
        if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,'\',strrep(shape_names{j},' ','_'),'_','plot','_N=', int2str(N),'.pdf')); 
        end
      
      fig_numbering=fig_numbering+1;
   end
   
   %plotting the first two coordinates of the samples
   if plot_sample_coordinates
       fig=figure(fig_numbering);
       for k=1:length(tags)
           samples=eval(strcat('samples_', tags{k}));
           %subplot(1,length(tags),k)
           root_ceil=ceil(sqrt(length(tags)));
           subplot(root_ceil,root_ceil,k);
           scatter(samples(1,:),samples(2,:),1,'Marker','.')
           title(title_tags{k});
           hold on; 
       end
       sgtitle(strcat('First Two Coordinates of Samples--',shape_names{j}))
       hold off;
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,'\',strrep(shape_names{j},' ','_'),'_','samples','_N=', int2str(N),'.pdf')); 
        end
       fig_numbering=fig_numbering+1;
   end
   
   %finding the r-values of the samples
   if plot_histograms|| rejection_rates
       dsts=zeros(max(N,N)*reps,k);
       for k=1:length(tags)
           if strcmp(tags{k},'adaptive')
               len=N;
           else
               len=N;
           end
           name=strcat('samples_',tags{k});
           mat=eval(name);
           for cc=1:reps
               for rr=1:len
                   i=(cc-1)*len+rr;
                   dsts(i,k)=dist(mat(1:d,i));
               end
           end
       end
   end
   
   %making a histogram of the r-values
   if plot_histograms
       fig=figure(fig_numbering);
       for k=1:length(tags)
           subplot(1,length(tags),k)
           if ~fixed_bins
               histogram(dsts(:,k),'Normalization','pdf');
           else
               histogram(dsts(:,k),num_bins,'Normalization','pdf');
           end
          
           title(title_tags{k});
           xlabel('r-values');
           ylabel('frequency');
       end
       sgtitle(strcat('histogram of radius values--', shape_names{j}))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,'\',strrep(shape_names{j},' ','_'),'_','r_hist','_N=', int2str(N),'.pdf')); 
       end
       fig_numbering=fig_numbering+1;
   end
   
   %making a histogram of the rejection rates
   if rejection_rates
       fig=figure(fig_numbering);
       disp(shape_names{j});
       for k=1:length(tags)
           name_acc=strcat('accepted_',tags{k});
           accepted=eval(name_acc);
           subplot(1,length(tags),k)
           if ~fixed_bins
               histogram(dsts(~accepted,k),'Normalization','count');
           else
               histogram(dsts(~accepted,k),num_bins,'Normalization','count');
           end
           title(title_tags{k});
           xlabel('r-values');
           ylabel('rejection rate');
           disp(strcat('acceptance rate, ',tags{k},': ', num2str(sum(accepted)/(N*reps))));
       end
       sgtitle(strcat('Rejection rates--', shape_names{j}))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gca, strcat(folder,'\',strrep(shape_names{j},' ','_'),'_','rejection_hist','_N=', int2str(N),'.pdf')); 
       end
       fig_numbering=fig_numbering+1;
       disp('--')
   end
   
end











function [id] =H(x,y,k)%the ratio of two proposal moves
    id=1;

end
function [y]=step_proposal(x,k,step_size)
    ell=length(k);
    delta=step_size*(2*rand(ell,1)-1);
    y=x;
    y(k)=y(k)+delta;
    y=mod(y,1);
end


function [k_next]=nxt(k_begin,dim)
    ell=length(k_begin);
    k_next=k_begin+ell-1;
    k_next=mod(k_next,dim)+1;
end





function[dst]=rect_then_square_cont_nonlinear(x,level)
    sz=size(x);
    d=sz(1);
    n=sz(2);
    c1=.5*ones(d,1);
    centers=repmat(c1,1,n);
    diff=x-centers;
    dist_sq=vecnorm(diff,inf);
    dist_sq_last=vecnorm(diff(2:d,:),inf);
    dist_first=abs(diff(1,:));
    k=1-level;
    dist_first=(-k+sqrt(k^2+4*dist_first))/2;
    dist_rect=max(dist_first,dist_sq_last);
    ind=dist_sq>=level;
    dst=zeros(length(ind),1);
    dst(ind)=dist_sq(ind);
    dst(~ind)=dist_rect(~ind);
    
    
end



function[dst]=rect_then_square_cont(x,level)
    sz=size(x);
    d=sz(1);
    n=sz(2);
    c1=.5*ones(d,1);
    centers=repmat(c1,1,n);
    diff=x-centers;
    dst_sq=vecnorm(diff,inf);
    dst_rect=vecnorm(diff(2:d,n),inf);
    ind=dst_sq>=level;
    dst=zeros(length(ind),1);
    dst(ind)=dst_sq(ind);
    dst(~ind)=dst_rect(~ind);
   
end


function[dst]=square(x)
    sz=size(x);
    d=sz(1);
    n=sz(2);
    c1=.5*ones(d,1);
    centers=repmat(c1,1,n);
    dst_sq=vecnorm(x-centers,inf);
    %dst_sq=0.5*ones(1,length(dst_sq))-dst_sq;
    dst=dst_sq;
   
end

function [rounded]=discretize(x,r)
    index=find(r>x,1,'first');
    rounded=r(index);
end