rng('shuffle')


ns=[4];%the number of disks
N=1000;%number of samples
reps=10;%number of repitions
torus=60;%the 60 or 90 degree torus
step_size=1;%the step size for the proposal
xtype='radius';%what to put on the x-axis: radius or density
yscale='linear';%log or linear for yscale
W=@(r) r; %the importance weight function. can be either a function or a numeric vector with num_points entries
num_points=20;%number of points in our discretization of r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters for plots and data saving
volume_vs_r=true; %do you want to see the plot of volume vs. r?
plot_sample_coordinates=true; %do you want to see the plot of the first two coordinates of the samples?
plot_r_histograms=true; %do you want to see a histogram of the rs?
rejection_frequencies=true; %do you want to see a histogram of rejection frequencies?
fig_numbering=1; %matlab starts numbering figures at this number
save_figures=true;%if we want to save the figures
save_data=true;%if we want to save the data
base_filename='run1';%base for all filenames
fixed_bins=true;%if we want a fixed number of bins in our histograms
num_bins=15;% the number of bins we want in our histograms


%parameters you want to change are probably above this line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%the following are dummy variables, should augment or remove
marginal=true;
methods_list=[marginal];
methods_names={'marginal'};
methods_number=sum(methods_list);
plot_labels={'marginal method'};
all_colors=getColorSet(1);
colors=all_colors(1,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return_samples=plot_sample_coordinates || plot_r_histograms || rejection_frequencies;
parent_directory='plots_and_data';%we put the directory 'folder' in this directory
folder='disks_partition_function';%location for saving figures and data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%here we check if the relevant folders exist and creat them if they don't
if ispc% file divider on windows vs unix and apple
    slash='\';
else
    slash='/';
end
if ~isfolder(parent_directory)
    mkdir(parent_directory);
end
folder=strcat(parent_directory,slash,folder);%TODO check that this and following 4 folder lines also runs on linux
if ~isfolder(folder)
    mkdir(folder);
end
data_folder=strcat(folder,slash,'data');
if save_data&&~isfolder(data_folder)
    mkdir(data_folder);
end
if save_data
    %saving the parameters which created the data
    save(strcat(data_folder,slash,base_filename,'parameters'));
end
if save_figures
    volumes_folder=strcat(folder,slash,'volume_plots');
    samples_folder=strcat(folder,slash,'samples_plots');
    r_histograms_folder=strcat(folder,slash,'r_histogram_plots');
    rejection_frequency_folder=strcat(folder,slash,'rejection_frequency_plots');
    if volume_vs_r && ~isfolder(volumes_folder)
        mkdir(volumes_folder);
    end
    if plot_sample_coordinates && ~isfolder(samples_folder)
        mkdir(samples_folder);
    end
    if plot_r_histograms && ~isfolder(r_histograms_folder)
        mkdir(r_histograms_folder);
    end
    if rejection_frequencies && ~isfolder(rejection_frequency_folder)
        mkdir(rejection_frequency_folder);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=length(ns)
    n=ns(j);
    %initializing the struct for storing the data
    for k=1:methods_number
        S(k).n=n;
        S(k).method=methods_names{k};
        S(k).plot_label=plot_labels{k};
        S(k).volumes=zeros(num_points,reps);%for holding the data
        if plot_sample_coordinates || plot_r_histograms || rejection_frequencies
            S(k).samples=zeros(2*n,reps*N);
        end
        if rejection_frequencies
            S(k).accepted=false(N*reps,1);
        end
        if plot_r_histograms|| rejection_frequencies
            S(k).samples_radiuses=zeros(N*reps,1);
        end
    end
    if torus==60%discretizing r
        [~,critical] = get_rigid_configurations(n,xtype,'tor60_on_90');
        r_max=max(critical);
    else
        r_max=sqrt(2)/4;
    end
    discretize_x_axis=linspace(0,r_max,num_points);%discretizing x axis
    if strcmp(xtype,'radius')
        r=discretize_x_axis;
        S(k).r=r;
    else
        if torus==90
            volume_torus=1;
        else
            volume_torus=(sqrt(3/2))^n;
        end
        S(k).density=discretize_x_axis;
        r=sqrt(discretize_x_axis/(n*pi));
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %generating the data
    
    for i=1:reps
        for k=1:methods_number
            method=S(k).method;
            if ~plot_sample_coordinates && ~plot_r_histograms && ~rejection_frequencies
                if strcmp(method,'marginal')
                    S(k).volumes(:,i)=disks_partition_function(N,n,r,torus,step_size,return_samples,W);
                end
            else
                if strcmp(method,'marginal')
                    [S(k).volumes(:,i),S(k).samples(:,((i-1)*N+1):(i*N)),~,~,S(k).accepted(((i-1)*N+1):(i*N))]=disks_partition_function(N,n,r,torus,step_size,return_samples,W);
                end
            end
        end
    end
   means=zeros(num_points,1);
   standard_deviations=zeros(num_points,1);
   for k=1:methods_number
      mat=S(k).volumes;
      means=mean(mat,2);
      standard_deviations=std(mat,1,2);
      S(k).means=means;
      S(k).standard_deviations=standard_deviations;
   end
   if plot_r_histograms|| rejection_frequencies
       for k=1:methods_number
            for i=1:reps*N
                sample=S(k).samples(:,i);
                if torus==90
                    dist=@dist_90;
                else
                    dist=@dist_60;
                end
                S(k).samples_radiuses(i)=dist(sample);
            end
       end
   end
   %saving the data
   if save_data
       save(strcat(data_folder,slash,base_filename,'n=',num2str(n)),'S');
   end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plotting
   %ploting volume vs. r
   if volume_vs_r
       fig=figure(fig_numbering);
       set(gca,'YScale', yscale);
       if strcmp(xtype,'radius')
           xx=r;
       else
           xx=density;
       end
       if torus==60
           len=length(critical)-1;
           if len>0
               handles_length=length(methods_names)+2;
           else
               handles_length=length(methods_names)+1;
           end
       else
           handles_length=length(methods_names);
       end
       plot_handles=gobjects(handles_length,1);
       legend_names=cell(handles_length,1);
       for k=1:length(methods_names)
           legend_names{k}=S(k).plot_label;
           means=S(k).means;
           standard_deviations=S(k).standard_deviations;
           plot_handles(k)=plot(xx,means,'Color', colors(k,:));
           hold on;
           plot(xx,means-2*standard_deviations,'Color',colors(k,:),'LineStyle', '--');
           hold on;
           plot(xx,means+2*standard_deviations,'Color',colors(k,:),'LineStyle', '--');
           hold on;
       end
       %plotting the critical densities/radiuses
       if ~strcmp(yscale,'log')
           mult=0;
       else
           mult=min(ylim);
       end
       if torus==60
           mx=max(critical);
           ind= critical==mx;
           plot_handles(length(methods_names)+1)=scatter(mx,mult, 50,"red",'*');
           hold on;
           v=critical(~ind);
           len=length(v);
           if len>0
               plot_handles(length(methods_names)+2)=scatter(v,mult*ones(1,len),50,"green",'*');
               hold on;
           end
       end
       %setting the legend
       if strcmp(xtype,'radius')
           loc='southwest';
       else
           if n<=4
               loc='southwest';
           else
               loc='northeast';
           end
       end
       
       if torus==60
           if strcmp(xtype, 'radius')
               legend_names{length(methods_names)+1}='maximum radius';
               if len>0
                   legend_names{length(methods_names)+2}='critical radiuses';
               end
           else
               legend_names(length(methods_names)+1)='maximum density';
               if len>0
                   legend_names(length(methods_names)+2)='critical densities';
               end
           end
       end
       legend(plot_handles,legend_names,'Location',loc);
       xlabel(xtype);
       ylabel('volume');
       title(strcat('volume vs. r--n=',num2str(n)));
       hold off;
       
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(volumes_folder,slash,base_filename,'_torus_',num2str(torus),'_','n=',num2str(n),'_N=', int2str(N),'.pdf'));
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
       sgtitle(strcat('First Two Coordinates of Samples--n=',num2str(n)))
       hold off;
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(samples_folder,slash,base_filename,'_torus_',num2str(torus),'_','n=',num2str(n),'_N=', int2str(N),'.pdf'));
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
       sgtitle(strcat('histogram of radius values--n=', num2str(n)))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(r_histograms_folder,slash,base_filename,'_torus_',num2str(torus),'_','n=',num2str(n),'_N=', int2str(N),'.pdf'));
       end
       fig_numbering=fig_numbering+1;
   end
   
   %making a histogram of the rejection frequencies
   if rejection_frequencies
       fig=figure(fig_numbering);
       disp(strcat('n=',num2str(n)));
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
           disp(strcat('acceptance rate, ',S(k).method,': ', num2str(sum(accepted)/(N*reps))));
       end
       sgtitle(strcat('Rejections--n=', num2str(n)))
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(rejection_frequency_folder,slash,base_filename,'_torus_',num2str(torus),'_','n=',num2str(n),'_N=', int2str(N),'.pdf'));
       end
       fig_numbering=fig_numbering+1;
       disp('--')
   end
    
end

