%rng('shuffle');
rng(5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the possible M(x)'s for which we want plots
square=false;%M(x)\leq r is an infinity ball radius r centered at .5 times the ones vector
rectangle_then_square=false;%for a parameter 'level' (defined later on),  M(x)\leq r is 
                            %1. if r\leq level, a hyperrectangle centered at 0.5 times the ones vector with one side length 2*level and the other side lengths 2*r
                            %2. if r\geq level, an infinity ball centered at 0.5 times the ones vector for r\geq level
rectangle_then_square_nonlinear=false;% %for a parameter 'level' (defined later on),  M(x)\leq r is
                            %1. if r\leq level, a hyperrectangle centered at 0.5 times the ones vector with one side length 2*(r*(1-level)+r^2) and the other side lengths 2*r. The side length 2*(r*(1-level)+r^2) was chosen to make the volume continuous but nonlinear on [0,level]
                            %2. if r\geq level, an infinity ball centered at 0.5 times the ones vector for r\geq level
disks90=false;%M(x)\leq r the configuration space for disks radius r in a 90 degree torus
disks60=true;%M(x)\leq r is the configuration space for disks radius r in a 60 degree torus 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters for plots
volume_vs_r=true; %do you want to see the plot of volume vs. r?
plot_sample_coordinates=true; %do you want to see the plot of the first two coordinates of the samples?
plot_r_histograms=true; %do you want to see a histogram of the rs?
rejection_frequencies=true; %do you want to see a histogram of rejection frequencies?
fig_numbering=4; %matlab starts numbering figures at this number
save_figures=true;%if we want to save the figures
fixed_bins=true;%if we want a fixed number of bins in our histograms
num_bins=15;% the number of bins we want in our histograms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the following variables toggle which methods we will be comparing
marginal_uniform_method=false; %uniform W for the marginal method
marginal_inverse_derivative_method=false; %weight function proportional to inverse derivative for the marginal method
marginal_antiderivative_inverse_volume_method=false; %weight function proportional to the antiderivative of the inverse weights for the marginal method
umbrella_antiderivative_inverse_volume_method=true; %weight function proportional to the antiderivative of the inverse weights for umbrella sampling

joint_antiderivative_inverse_volume_method=false;%if we want to include the method that samples in the joint distribution
joint_inverse_derivative_method=false;
adaptive_antiderivative_inverse_volume_method=false;%adaptive method converging to weights proportional to the antiderivative of the inverse weights
wang_landau_method=false;%if we want to include the wang-landau method



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters for methods-- one we may want to change
d=4;%the dimension
num_points=51;%number of points in our discretization
N=100000;%number of samples
burn_in=100;%burn_in_time
reps=10;%number of repetitions
step_size=1; %note--optimal for wang-landau in 2d seems to be .1
level=.25;%for the rectangle_then_square/rectangle_then_square_nonlinear shapes, the r at which M(x)\leq r switches from being a rectangle to a square
forget_rate_anti_adaptive=10^-8;%rate parameter for the inverse volumes method
p=0.1;%rate of changing r in joint distribution.



%parameters you want to change are probably above this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parameters for methods that we probably don't want to change around
proposal=@(x,k) step_proposal(x,k,step_size);%the proposal function

%for proposing which coordinates to move
next=@(k)nxt(k,d);%picking the next coordinates to move in the proposal
start_coordinates=[1];%the first coordinates to move in the proposal. for the umbrella sampling methods, only 1-dimensional 'start_coordinates' has been implemented



%some extra parameters for the Wang-Landau method
WL_dist_discrete=false;%if we want to use the discrete version of the WL method
log_scale_calculations=true;%log scale calculations in the wang-landau method
tolerance=.2; %the some of the wang landau updates use histograms of M-values. This tolerance controls how close the histogram should be to uniform before updating f 
f_update=@(f,g,N,Hist,r,dst,other)f_update_optimal(f,g,N,Hist,r,dst,other,log_scale_calculations,tolerance);
%f_update=@(f,g,N,Hist,r,dst,log_scale,other)f_update_inverse_time(f,g,N,Hist,r,dst,log_scale,other);


%weight parameters for our methods
W_uniform=ones(num_points,1);%uniform weights
w0=ones(num_points,1); %starting point for adaptive weights

folder='simple_volumes_figures';
parent_directory='plots_and_data';%we put the directory 'folder' in this directory
if ispc% file divider on windows vs unix and apple
    slash='\';
else
    slash='/';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%listing the M-functions which determine our shapes and the shape names

%handles for our shapes
square_handle=@(x)square_function(x);
rect_square_handle=@(x)rect_then_square(x,level);
rect_square_nonlinear_handle=@(x)rect_then_square_nonlinear(x,level);

M_function_list={};
shape_names={};
if square
    len=length(shape_names);
    M_function_list{len+1}=square_handle;
    shape_names{len+1}='square';
end
if rectangle_then_square
    len=length(shape_names);
    M_function_list{len+1}=rect_square_handle;
    shape_names{len+1}='rectangle then square';
end
if rectangle_then_square_nonlinear
    len=length(shape_names);
    M_function_list{len+1}=rect_square_nonlinear_handle;
    shape_names{len+1}='rectangle then square nonlinear';
end
if disks90
    len=length(shape_names);
    M_function_list{len+1}=@dist_90;
    shape_names{len+1}='disks90';
end
if disks60
    len=length(shape_names);
    M_function_list{len+1}=@dist_60;
    shape_names{len+1}='disks60';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%labeling our methods and assigning colors
all_methods_booleans=[marginal_uniform_method; marginal_inverse_derivative_method; marginal_antiderivative_inverse_volume_method; adaptive_antiderivative_inverse_volume_method; wang_landau_method; joint_antiderivative_inverse_volume_method; joint_inverse_derivative_method; umbrella_antiderivative_inverse_volume_method; ];
methods_list={};
plot_labels={};
methods_number=sum(all_methods_booleans);%number of stochastic methods
all_colors=getColorSet(length(all_methods_booleans)+1);%+1 for analytically finding the volumes
colors=all_colors(1,:);
color_index=2;
if marginal_uniform_method
   len=length(methods_list);
   methods_list{len+1}='uniform';
   plot_labels{len+1}='uniform weights';
   colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if marginal_inverse_derivative_method
   len=length(methods_list);
   methods_list{len+1}='inverse_derivative';
   plot_labels{len+1}='inverse derivative weights';
   colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if marginal_antiderivative_inverse_volume_method
   len=length(methods_list);
   methods_list{len+1}='anti_inverse_volume';
   plot_labels{len+1}='integral inverse volume weights';
   colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if adaptive_antiderivative_inverse_volume_method
   len=length(methods_list);
   methods_list{len+1}='adaptive_anti_inverse_volume';
   plot_labels{len+1}='adaptive integral inverse volume';
   colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if wang_landau_method
    len=length(methods_list);
    methods_list{len+1}='wang_landau';
    plot_labels{len+1}='Wang-Landau method';
    colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if joint_antiderivative_inverse_volume_method
    len=length(methods_list);
    methods_list{len+1}='joint_anti_inverse_volume';
    plot_labels{len+1}='joint integral inverse volume';
    colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if joint_inverse_derivative_method
    len=length(methods_list);
    methods_list{len+1}='joint_inverse_derivative';
    plot_labels{len+1}='joint inverse derivative';
    colors(len+2,:)=all_colors(color_index,:);
end
color_index=color_index+1;
if umbrella_antiderivative_inverse_volume_method
    len=length(methods_list);
    methods_list{len+1}='umbrella_anti_inverse_volume';
    plot_labels{len+1}='umbrella integral inverse volume';
    colors(len+2,:)=all_colors(color_index,:);
end
    
%initializing the struct for storing the data
for k=1:length(methods_list)
    %labeling the methods
    method=methods_list{k};
    S(k).method=method;
    S(k).plot_label=plot_labels{k};
    %initializing matricies that will store volumes
    S(k).volumes=zeros(num_points,reps);
    
    %initializing matrices to store samples
    %and vectors that will store booleans representing if a proposal
    %was accepted or not
    if plot_sample_coordinates || plot_r_histograms || rejection_frequencies
        if strcmp(method, 'joint_anti_inverse_volume')|| strcmp(method, 'joint_inverse_derivative')
            S(k).samples=zeros(d+1,reps*N);
        else
            S(k).samples=zeros(d,reps*N);
        end
        S(k).accepted=false(N*reps,1);
    end


    
    if plot_r_histograms|| rejection_frequencies
        S(k).samples_radiuses=zeros(N*reps,1);
    end
    
end


if ~isfolder(parent_directory)&&save_figures
    mkdir(parent_directory);
end
folder=strcat(parent_directory,slash,folder);%TODO check that this line also runs on unix
if ~isfolder(folder)&&save_figures
    mkdir(folder);
end

for j=1:length(M_function_list)
   dist=M_function_list{j};
   exact_volumes=zeros(num_points,1);
   volume_derivative=zeros(num_points,1);
   %variables unique for each shape
   if isequal(dist,square_handle)
      reference_volume=1;%reference volume
      r=linspace(0,0.5,num_points)'; %discretizing radiuses
      up=true;%our reference volume is the largest shape
      increasing=true;
      exact_volumes=(2*r).^d;
      anti_inverse_volume_weights=r.^(-d+1);
      volume_derivative=2*d*(2*r).^(d-1);
      r_lim=[0,1/2];%upper and lower limits on r
      proposal_umbrella=@square_hit_and_run;
   elseif isequal(dist,rect_square_handle)
       reference_volume=1;%reference volume
       r=linspace(0,0.5,num_points)'; %discretizing radiuses
       up=true;%our reference volume is the largest volume
       increasing=true;
       ind=r<=(level);
       %calculating the exact volumes and volume derivatives
       exact_volumes(ind)=(2*level).*(2*(r(ind))).^(d-1);
       volume_derivative(ind)=(2*level)*2*(d-1)*(2*r(ind)).^(d-2);
       exact_volumes(~ind)=(2*(r(~ind))).^d;
       volume_derivative(~ind)=2*d*(2*r(~ind)).^(d-1);
       r_lim=[0,1/2];%upper and lower limits on r
       proposal_umbrella=@(x,k,m)rect_then_square_hit_and_run(x,k,m,level);
   elseif isequal(dist,rect_square_nonlinear_handle)
       reference_volume=1;%reference volume
       r=linspace(0,0.5,num_points)'; %discretizing radiuses
       up=true;%our reference volume is the largest volume
       increasing=true;
       %calculating the exact volumes and volume derivatives
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
       r_lim=[0,1/2];%upper and lower limits on r
       proposal_umbrella=@(x,k,m)rect_then_square_nonlinear_hit_and_run(x,k,m,level);
   elseif isequal(dist, @dist_90)
       reference_volume=1;%reference volume
       r_max=sqrt(2)/4;%discretizing the radiuses
       r=linspace(0,r_max,num_points)';
       if d~=4%the exact volumes are implemented only for 2 disks
           error('for disks in a torus, exact volumes only known for 2 disks')
       end
       up=true;%our reference volume is the largest volume
       increasing=false;
       %calculating the exact volumes and volume derivatives
       ind=r>1/4;
       exact_volumes(~ind)=1-4*pi*r(~ind).^2;
       exact_volumes(ind)=1-4*pi*r(ind).^2+16*r(ind).^2 .*acos(1./(4*r(ind)))-sqrt(16*r(ind).^2-1);
       volume_derivative= 8*pi*r;
       r_above=r(ind);
       volume_derivative(ind)= volume_derivative(ind)-32*r_above.*acos(1./(4*r_above))-4./sqrt(1-1./(4*r_above).^2)+16*r_above./sqrt(16*r_above.^2-1);
       volume_derivative=max(volume_derivative,0);
       r_lim=[0,r_max];%upper and lower limits on r
       A=eye(2,2);
       proposal_umbrella=@(x,k,m)disks_single_coordinate_hit_and_run(x,k,m,A);
   elseif isequal(dist,@dist_60)
       reference_volume=sqrt(3)/2;%reference volume
       if d~=4
           error('for disks in a torus, exact volumes only known for 2 disks')
       end
       up=true;%the exact volumes are implemented only for 2 disks
       increasing=false;
       [xs,radiuses] = get_rigid_configurations(d/2,'radius','tor60_on_90');
       r_max=max(radiuses);
       r=linspace(0,r_max,num_points)';%discretizing radiuses
       %calculating the exact volumes and volume derivatives
       ind=r>1/4;
       exact_volumes(~ind)=sqrt(3)/2-4*pi*r(~ind).^2;
       exact_volumes(ind)=sqrt(3)/2-4*pi*r(ind).^2+3/2*(16*r(ind).^2 .*acos(1./(4*r(ind)))-sqrt(16*r(ind).^2-1));
       volume_derivative= 8*pi*r;
       r_above=r(ind);
       volume_derivative(ind)=volume_derivative(ind)+3/2*(-32*r_above.*acos(1./(4*r_above))-4./sqrt(1-1./(4*r_above).^2)+16*r_above./sqrt(16*r_above.^2-1));
       volume_derivative(1)=0;
       volume_derivative(num_points)=0;
       r_lim=[0,r_max];%limits on r
       A=[1 1/2;
          0 sqrt(3)/2];
       proposal_umbrella=@(x,k,m)disks_single_coordinate_hit_and_run(x,k,m,A);
   end
   if WL_dist_discrete
        wdist=@(x)discretize(dist(x),r);
   else
        wdist=dist;
   end
   inverse_derivative=1./volume_derivative;
   quadrature_weights=non_uniform_trapezoidal_weights(r);
   inverse_volumes=1./exact_volumes;
   indd=inverse_volumes==inf;
   inverse_volumes(indd)=max(inverse_volumes(~indd));
   if increasing
      anti_inverse_volume_weights=flip(cumsum(quadrature_weights./flip(inverse_volumes)));
   else
      anti_inverse_volume_weights=cumsum(quadrature_weights./inverse_volumes); 
   end
   Ns_anti_inverse_volumes=get_Ns_umbrella(N,r,increasing, anti_inverse_volume_weights);
   burn_ins_anti_inverse_volumes=zeros(length(Ns_anti_inverse_volumes),1);
   if increasing
        burn_ins_anti_inverse_volumes(1)=burn_in;
   else
       burn_ins_anti_inverse_volumes(length(Ns_anti_inverse_volumes))=burn_in;
   end
   if umbrella_antiderivative_inverse_volume_method
       if isequal(dist,@dist_60)||isequal(dist,@dist_90)
           dist_umbrella=@(x)dist(x(1:(length(x)-2)));
           H_umbrella=@(y,x,k,m)H_umbrella_disks(y,x,k,m,A);
       else
           dist_umbrella=dist;
           H_umbrella=@(y,x,k,m)1;
       end
   end
   %looping over methods
   for i=1:reps
        x0=rand(d,1);
        if isequal(dist,@dist_60)||isequal(dist,@dist_90)
            x0_umbrella=[x0;0;0];
        else
            x0_umbrella=x0;
        end
        for k=1:length(methods_list)
            method=S(k).method;
            if ~plot_sample_coordinates && ~plot_r_histograms && ~rejection_frequencies
                if strcmp(method,'uniform')
                    S(k).volumes(:,i)=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,W_uniform,r);
                end
                if strcmp(method,'anti_inverse_volume')
                    S(k).volumes(:,i)=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,anti_inverse_volume_weights,r);
                end
                if  strcmp(method,'inverse_derivative')
                    S(k).volumes(:,i)=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,inverse_derivative,r);
                end
                if strcmp(method,'adaptive_anti_inverse_volume')
                    S(k).volumes(:,i)=volume_adaptive_integral_inverse_volumes(N,burn_in,increasing,x0,next, start_coordinates,  proposal, @H,r,dist,plot_sample_coordinates,reference_volume,forget_rate_anti_adaptive, w0);
                end
                if strcmp(method,'wang_landau')
                    S(k).volumes(:,i)=volume_wang_landau(N,increasing,x0,start_coordinates,next,proposal,@H,f_update,r,wdist,reference_volume, plot_sample_coordinates,log_scale_calculations);
                end
                if strcmp(method,'joint_anti_inverse_volume')
                    S(k).volumes(:,i)=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,anti_inverse_volume_weights,[]);
                end
                if strcmp(method,'joint_inverse_derivative')
                    S(k).volumes(:,i)=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist, plot_sample_coordinates,reference_volume,inverse_derivative,[]);
                end
                if strcmp(method,'umbrella_anti_inverse_volume')
                   S(k).volumes(:,i)=volume_umbrella_sampling(Ns_anti_inverse_volumes,burn_ins_anti_inverse_volumes,up,increasing,x0_umbrella, next,start_coordinates, proposal_umbrella,H_umbrella,r,dist_umbrella, false,reference_volume); 
                end
            else
                if strcmp(method,'uniform')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~,S(k).accepted(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing, x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,W_uniform,r);
                end
                if strcmp(method,'anti_inverse_volume')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~,S(k).accepted(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,anti_inverse_volume_weights,r);
                end
                if strcmp(method,'inverse_derivative')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~,S(k).accepted(((i-1)*N+1):(i*N))]=volume_marginal(N,burn_in,up,increasing,x0, next,start_coordinates, proposal,@H,r,dist, true,reference_volume,inverse_derivative,r);
                end
                if strcmp(method,'adaptive_anti_inverse_volume')
                    [S(k).volumes(:,i),~,~,S(k).samples(:,((i-1)*N+1):(i*N)),S(k).accepted(((i-1)*N+1):(i*N))]=volume_adaptive_integral_inverse_volumes(N,burn_in,increasing, x0,next, start_coordinates,  proposal, @H,r,dist,true,reference_volume,forget_rate_anti_adaptive, w0);
                end
                if strcmp(method,'wang_landau')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)), S(k).accepted(((i-1)*N+1):(i*N)),~]=volume_wang_landau(N,increasing,x0,start_coordinates,next,proposal,@H,f_update,r,wdist,reference_volume,true,log_scale_calculations);
                end
                if strcmp(method,'joint_anti_inverse_volume')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~, S(k).accepted(((i-1)*N+1):(i*N))]=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist,true,reference_volume,anti_inverse_volume_weights,[]);
                end
                if strcmp(method,'joint_inverse_derivative')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~, S(k).accepted(((i-1)*N+1):(i*N))]=volume_joint(N,burn_in,up,increasing,x0, next,start_coordinates,p, proposal,@H,r,dist,true,reference_volume,inverse_derivative,[]);
                end
                if strcmp(method,'umbrella_anti_inverse_volume')
                   [S(k).volumes(:,i),smp,~,acpt]=volume_umbrella_sampling(Ns_anti_inverse_volumes,burn_ins_anti_inverse_volumes,up,increasing,x0_umbrella, next,start_coordinates, proposal_umbrella,H_umbrella,r,dist_umbrella, true,reference_volume); 
                   accepted_temp=false(N,1);
                   if isequal(dist,@dist_60)||isequal(dist,@dist_90)                      
                       samples_temp=zeros(d+2,N);
                   else
                       samples_temp=zeros(d,N);
                   end
                   bottom_index=1;
                   for ell=1:length(Ns_anti_inverse_volumes)
                        N_temp=Ns_anti_inverse_volumes(ell); 
                        top_index=bottom_index+N_temp-1;
                        samples_temp(:,bottom_index:top_index)=smp{ell};
                        accepted_temp(bottom_index:top_index)=acpt{ell};
                        bottom_index=top_index+1;
                   end
                   S(k).accepted(((i-1)*N+1):(i*N))=accepted_temp;
                   S(k).samples(:,((i-1)*N+1):(i*N))=samples_temp(1:d,:);
                end
            end
        end
   end
    
   %averaging
   means=zeros(num_points,1);
   standard_deviations=zeros(num_points,1);
   for k=1:methods_number
      mat=S(k).volumes;
      means=mean(mat,2);
      standard_deviations=std(mat,1,2);
      S(k).means=means;
      S(k).standard_deviations=standard_deviations;
   end
   %calculating the maximum radius corresponding to each sample
   if plot_r_histograms|| rejection_frequencies
       for k=1:methods_number
            for i=1:reps*N
                sample=S(k).samples(:,i);
                method=S(k).method;
                if strcmp(method,'joint_anti_inverse_volume')|| strcmp(method,'joint_inverse_derivative')
                    sample=sample(1:d);
                end
                S(k).samples_radiuses(i)=dist(sample);
            end
       end
   end
   
   %ploting volume vs. r
   if volume_vs_r
      fig=figure(fig_numbering);
      plot_handles=gobjects(methods_number+1,1);
      legend_names=cell(methods_number+1,1);
      plot_handles(1)=plot(r, exact_volumes,'Color',colors(1,:));
      hold on;
      legend_names{1}='true volume';
      for k=2:length(methods_list)+1
         means=S(k-1).means;
         standard_deviations=S(k-1).standard_deviations;
         legend_names{k}=S(k-1).plot_label;
         plot_handles(k)=plot(r,means,'Color', colors(k,:));
         hold on;
         plot(r,means-2*standard_deviations,'Color',colors(k,:),'LineStyle', '--');
         hold on;
         plot(r,means+2*standard_deviations,'Color',colors(k,:),'LineStyle', '--');
         hold on;
      end
        if increasing
            loc='Southeast';
        else
            loc='Southwest';
        end
        legend(plot_handles,legend_names,'Location',loc);
        xlabel('r');
        ylabel('volume');
        title(strcat('volume vs. r--',shape_names{j}));
        hold off;
          
        if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,slash,strrep(shape_names{j},' ','_'),'_','plot','_N=', int2str(N),'.pdf')); 
        end
      
      fig_numbering=fig_numbering+1;
   end
   
   %plotting the first two coordinates of the samples
   if plot_sample_coordinates
       fig=figure(fig_numbering);
       
       for k=1:length(methods_list)
           samples=S(k).samples;
           %subplot(1,length(tags),k)
           root_ceil=ceil(sqrt(length(methods_list)));
           subplot(root_ceil,root_ceil,k);
           scatter(samples(1,:),samples(2,:),1,'Marker','.','MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410])
           title(S(k).plot_label);
           hold on; 
       end
       sgtitle(strcat('First Two Coordinates of Samples--',shape_names{j}))
       hold off;
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,slash,strrep(shape_names{j},' ','_'),'_','samples','_N=', int2str(N),'.pdf')); 
        end
       fig_numbering=fig_numbering+1;
   end
   
   
   %making a histogram of the r-values
   if plot_r_histograms
       fig=figure(fig_numbering);
       for k=1:length(methods_list)
           subplot(1,length(methods_list),k)
           m_list=S(k).samples_radiuses;
           if ~fixed_bins
               hst=histogram(m_list,'Normalization','pdf');
           else
               hst=histogram(m_list,num_bins,'Normalization','pdf');
           end
           hst.FaceColor=[0, 0.4470, 0.7410];
          
           title(S(k).plot_label);
           xlabel('r-values');
           ylabel('frequency');
       end
       sgtitle(strcat('histogram of radius values--', shape_names{j}))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(folder,slash,strrep(shape_names{j},' ','_'),'_','r_hist','_N=', int2str(N),'.pdf')); 
       end
       fig_numbering=fig_numbering+1;
   end
   
   %making a histogram of the rejection rates
   if rejection_frequencies
       fig=figure(fig_numbering);
       disp(shape_names{j});
       for k=1:length(methods_list)
           accepted=S(k).accepted;
           m_list=S(k).samples_radiuses;
           subplot(1,length(methods_list),k)
           if ~fixed_bins
               hst=histogram(m_list(~accepted),'Normalization','count');
           else
               hst=histogram(m_list(~accepted),num_bins,'Normalization','count');
           end
           hst.FaceColor=[0, 0.4470, 0.7410];
           title(S(k).plot_label);
           xlabel('r-values');
           ylabel('rejection frequency');
           disp(strcat('acceptance rate, ',methods_list{k},': ', num2str(sum(accepted)/(N*reps))));
       end
       sgtitle(strcat('Rejections--', shape_names{j}))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gca, strcat(folder,slash,strrep(shape_names{j},' ','_'),'_','rejection_frequency_hist','_N=', int2str(N),'.pdf')); 
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


function [next_coordinates]=nxt(coordinates,dim)
    ell=length(coordinates);
    next_coordinates=coordinates+ell-1;
    next_coordinates=mod(next_coordinates,dim)+1;
end





function[dst]=rect_then_square_nonlinear(x,level)
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



function[dst]=rect_then_square(x,level)
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


function[dst]=square_function(x)
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
function[y]=disks_single_coordinate_hit_and_run(x,k,m,A)
    if k==1
        disk=[1;2];
        theta=0;
    elseif k==2
        disk=[1;2];
        theta=pi/2;
    elseif k==3
        disk=[3;4];
        theta=0;
    else
        disk=[3;4];
        theta=pi/2;
    end
%     d=length(x);
%     y=zeros(d,1);
%     x_disks=x(1:(d-2));
%     x_move=maximum_move_disks(x_disks,disk,m,theta,A); %calculating the maximum distance disk k could travel before hitting another disk, (bounded above by some constant)
%     
%     rr=x_move*(rand(1,1)); %randomly picking a distance to move
%     delta=rr*[cos(theta);sin(theta)];
%     y(1:(d-2))=x_disks;
%     y(disk)=mod(y(disk)+delta,1);%moving disk 'disk'
%     y(d-1)=theta;
%     y(d)=x_move; 
    y=components_proposal(x,disk,m,A);
end

function [ratio]=H_umbrella_disks(y,x,k,m,A)
    if k==1||k==2
        disk=[1 ;2];
    else
        disk=[3;4];
    end
    ratio=H_components_proposal(y,x,disk,m,A);
end

function[y]=rect_then_square_nonlinear_hit_and_run(x,k,m,level)
                            %1. if r\leq level, a hyperrectangle centered at 0.5 times the ones vector with one side length 2*(r*(1-level)+r^2) and the other side lengths 2*r. The side length 2*(r*(1-level)+r^2) was chosen to make the volume continuous but nonlinear on [0,level]
                            %2. if r\geq level, an infinity ball centered at 0.5 times the ones vector for r\geq level
    if m>=level
        y=square_hit_and_run(x,k,m);
    else
        d=length(x);
        y=x;
        if k~=1
            y(k)=2*m*(rand(1,1))+(.5-m);
        else
            side_length=2*(m*(1-level)+m^2);
            y(1)=side_length*rand(1,1)+.5-side_length/2;
            y(2:d)=x(2:d);
        end
    end
end
function [y]=rect_then_square_hit_and_run(x,k,m,level)
    if m>=level
       y=square_hit_and_run(x,k,m);
    else
        d=length(x);
        y=x;
        if k~=1
            y(k)=2*m*(rand(1,1))+(.5-m);
        else
            y(1)=2*level*rand(1,1)+.5-level;
            y(2:d)=x(2:d);
        end
    end
end

function [y]=square_hit_and_run(x,k,m)
    y=x;
    y(k)=2*m*(rand(1,1))+(.5-m);
end