rng('shuffle')


n=[4];%the number of disks
component_numbers='all';%can be a list of component numbers or 'all'
N=100;%number of samples
burn_in=100;
reps=10;%number of repitions
xtype='radius';%what to put on the x-axis: radius or density
yscale='linear';%log or linear for yscale
W=@(r) r; %the importance weight function. can be either a function or a numeric vector with num_points entries
num_points=20;%number of points in our discretization of r
upper_fraction=1;%can be any number in [0,1] we consider r in [upper_fraction*c,c], where c is the radius of the densest critical configuration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parameters for plots and data saving
volume_vs_r=true; %do you want to see the plot of volume vs. r?
diagnostic_plot_indices=5;%the radiuses at which we want diagnostic plots. Can be either 1. a list of integers between 1 and num_points, with 'auto_space_diagnostic_radiuses' as false
auto_space_diagnostic_radiuses=true;                                                   % 2. a single integer between 1 and num_points, with 'auto_space_diagnostic_radiuses' as true. diagnostic plots will then be created for 'diagnostic_plot_radiuses' equally spaced points                                                 
plot_sample_coordinates=true; %do you want to see the plots of the first two coordinates of the samples?
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
%%%%%%%%%%%%%%%%%%%%%%%%
return_samples=plot_sample_coordinates || plot_r_histograms || rejection_frequencies;
folder='disks_components_volumes';%location for saving figures and data
parent_directory='plots_and_data';%we put the directory 'folder' in this directory

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
    rejection_frequency_folder=strcat(folder,slash,'rejection_rate_plots');
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
%we now list out the components, and the corresponding initial values and
%radiuses
[critical_configurations,critical_radiuses] = get_rigid_configurations(n,xtype,'tor60_on_90');
if strcmp(component_numbers,'all')
    component_numbers=length(critical_radiuses);
elseif any(length(critical_radiuses)<component_numbers)
    error(strcat("there aren't that many components for n=",num2str(n), "disks"))
else
     critical_configurations=critical_configurations(:,component_numbers);
     critical_radiuses=critical_radiuses(component_numbers);
end

%determining the r's for which we make the diagnostic plots
if auto_space_diagnostic_radiuses
    if length(diagnostic_plot_indices)~=1 
        error("if 'auto_space_diagnostic_radiuses' is true, then  'diagnostic_plot_radiuses' should be a single integer")
    end
    if length(diagnostic_plot_indices)>num_points
       error("the number of radiuses for which diagnostic plots are created (diagnostic_plot_radiuses) should be at most the number of points in the discretization (num_points)")
    end
    if diagnostic_plot_indices==num_points
        diagnostic_plot_indices=1:num_points;
    elseif diagnostic_plot_indices==num_points-1
        diagnostic_plot_indices=1:(num_points-1);
    else
       equal_space_points=linspace(1,num_points,diagnostic_plot_indices+2);
       equal_space_points=equal_space_points(2:(diagnostic_plot_indices+1));
       diagnostic_plot_indices=round(equal_space_points);
    end
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializing the struct the holds our data
for j=length(component_numbers)
    component_number=component_numbers(j);
    critical_r=critical_radiuses(j);
    critical_configuration=critical_configurations(:,j);
    %initializing the struct for storing the data
    for k=1:methods_number
        S(k).n=n;
        S(k).component_number=component_number;
        S(k).method=methods_names{k};
        S(k).plot_label=plot_labels{k};
        S(k).volumes=zeros(num_points,reps);%for holding the data
        if plot_sample_coordinates || plot_r_histograms || rejection_frequencies
            S(k).samples=cell(num_points,1);
        end
        if rejection_frequencies
            S(k).accepted=false(num_points,N*reps);
        end
        if plot_r_histograms|| rejection_frequencies
            S(k).samples_radiuses=zeros(num_points,N*reps);
        end
    end
    discretize_x_axis=linspace((1-upper_fraction)*critical_r,critical_r,num_points);%discretizing x axis
    if strcmp(xtype,'radius')
        r=discretize_x_axis;
        S(k).r=r;
    else
        volume_torus=(sqrt(3/2))^n;
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
                    S(k).volumes(:,i)=volume_disks_components(N,burn_in,critical_configuration,r,60,return_samples, W);
                end
            else
                if strcmp(method,'marginal')
                    [S(k).volumes(:,i),samples,~,~,S(k).accepted(:,((i-1)*N+1):(i*N))]=volume_disks_components(N,burn_in,critical_configuration,r,60,return_samples, W);
                    for ell=1:num_points
                        if ell<num_points
                            S(k).samples{ell}(:,((i-1)*N+1):(i*N))=samples{ell};
                        else
                            S(k).samples{ell}(:,((i-1)*N+1):(i*N))=repmat(critical_configuration,1,N);
                        end
                    end
                    
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
           for ell=1:num_points
               samples=S(k).samples{ell};
               for i=1:reps*N
                   sample=samples(:,i);
                   S(k).samples_radiuses(ell,i)=dist_60(sample);
               end
           end
       end
   end
   %saving the data
   if save_data
       save(strcat(data_folder,slash,base_filename,'n=',num2str(n),'component=',num2str(component_number)),'S');
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
       handles_length=methods_number+1;    
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
       plot_handles(length(methods_names)+1)=scatter(critical_r,mult, 50,"red",'*');


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
       legend_names{length(methods_names)+1}='maximum radius';
       legend(plot_handles,legend_names,'Location',loc);
       xlabel(xtype);
       ylabel('volume');
       over_title=strcat('volume vs. r--n=',num2str(n),{', '}, 'component= ',num2str(component_number));
       title(over_title{1});
       hold off;
       
       if save_figures
           set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
           exportgraphics(gcf, strcat(volumes_folder,slash,base_filename,'_','n=',num2str(n),'component=',num2str(component_number),'_N=', int2str(N),'.pdf'));
       end
       
       fig_numbering=fig_numbering+1;
   end
   %plotting the first two coordinates of the samples
   if plot_sample_coordinates
       for k=1:length(methods_list)
           method_samples=S(k).samples;
           root_ceil=ceil(sqrt(length(diagnostic_plot_indices)));
           fig=figure(fig_numbering);
           for i=1:length(diagnostic_plot_indices)
               samples=method_samples{diagnostic_plot_indices(i)};
               r_current=r(diagnostic_plot_indices(i));
               subplot(root_ceil,root_ceil,i);
               scatter(samples(1,:),samples(2,:),1,'Marker','.','MarkerFaceColor',[0, 0.4470, 0.7410],'MarkerEdgeColor',[0, 0.4470, 0.7410])
               if strcmp(xtype, 'radius')
                   title(strcat('radius=',num2str(round(r_current,3))));
               else
                   title(strcat('density=',num2str(round(r_current,3))));
               end
               hold on;
               scatter(critical_configuration(1),critical_configuration(2),20,'Marker','.','MarkerFaceColor',[1, 0, 0],'MarkerEdgeColor',[1, 0, 0])
               hold on;
               num_discretization=500;
               t=linspace(0,2*pi,num_discretization);
               R=max(dist_60(critical_configuration)-r_current,0);
               x=mod(R*cos(t)+critical_configuration(1),1);
               y=mod(R*sin(t)+critical_configuration(2),1);
               scatter(x,y,1,'Marker','.','MarkerFaceColor',[1, 0, 0],'MarkerEdgeColor',[1, 0, 0]);
               hold on;
           end
           over_title=strcat(S(k).plot_label,', First Two Coordinates of Samples--n=',num2str(n),{', '},'component=',num2str(component_number));
           sgtitle(over_title);
           hold off;
           if save_figures
               set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
               exportgraphics(gcf, strcat(samples_folder,slash,base_filename,'_','n=',num2str(n),'component=',num2str(component_number),'_N=', int2str(N),'.pdf'));
           end
           fig_numbering=fig_numbering+1;
       end
   end
   
   
   %making a histogram of the r-values
   if plot_r_histograms 
       for k=1:length(methods_list)
           fig=figure(fig_numbering);
           r_lists=S(k).samples_radiuses;
           for i=1:length(diagnostic_plot_indices)
               r_list=r_lists(diagnostic_plot_indices(i),:);
               r_current=r(diagnostic_plot_indices(i));
               subplot(1,length(diagnostic_plot_indices),i)
               if ~fixed_bins
                   hst=histogram(r_list,'Normalization','pdf');
               else
                   hst=histogram(r_list,num_bins,'Normalization','pdf');
               end
              
               hst.FaceColor=[0, 0.4470, 0.7410];

               title(strcat('radius=',num2str(round(r_current,3))));
               ylabel('frequency');
           end
           over_title=strcat(S(k).plot_label,', histogram of radius values--n=', num2str(n),{', '},'component=',num2str(component_number));
           sgtitle(over_title{1})
           
           if save_figures
               set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
               exportgraphics(gcf, strcat(r_histograms_folder,slash,base_filename,'_','n=',num2str(n),'component=',num2str(component_number),'_N=', int2str(N),'.pdf'));
           end
           fig_numbering=fig_numbering+1;
           
       end
   end
   
   %making a histogram of the rejection rates
   if rejection_frequencies
       for k=1:length(methods_list)
           fig=figure(fig_numbering);
           disp(strcat(S(k).method,', component=',num2str(component_number)));
           accepted_list=S(k).accepted;
           r_lists=S(k).samples_radiuses;
           for i=1:length(diagnostic_plot_indices)
               r_list=r_lists(diagnostic_plot_indices(i),:);
               r_current=r(diagnostic_plot_indices(i));
               accepted=accepted_list(diagnostic_plot_indices(i),:);
               subplot(1,length(diagnostic_plot_indices),i)
               if ~fixed_bins
                   hst=histogram(r_list(~accepted),'Normalization','count');
               else
                   hst=histogram(r_list(~accepted),num_bins,'Normalization','count');
               end
               hst.FaceColor=[0, 0.4470, 0.7410];
               title(strcat('radius=',num2str(round(r_current,3))));
               xlabel('r-values');
               ylabel('rejection frequency');
               disp(strcat('acceptance rate, ',S(k).method,'component number:',component_number,' r',num2str(r_current),': ', num2str(sum(accepted)/(N*reps))));
           end
           over_title=strcat(S(k).plot_label,', Rejection frequency--n=', num2str(n),{' '},'component=',num2str(component_number));
           sgtitle(over_title{1})
           if save_figures
               set(fig, 'units', 'inches', 'position', [2 3 7 7])% [left bottom width height]
               exportgraphics(gcf, strcat(rejection_frequency_folder,slash,base_filename,'_','n=',num2str(n),'component=',num2str(component_number),'_N=', int2str(N),'.pdf'));
           end
           fig_numbering=fig_numbering+1;
           disp('--')
       end
   end
    
end






