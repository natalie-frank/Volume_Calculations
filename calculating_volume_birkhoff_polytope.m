%for calculating the volume of the birkhoff polytope
rng(10)
d_max=12;
N=1000;
burn_in=1000;
reps=10;
subspace_normalization=true;
return_samples=false;
node_num=[];

W=@(m) exp(-m);
%W=@(r)ones(length(r),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d_list=2:d_max;
birkhoff_exact_volumes=[ log(2), log(9/8), log(176)-log(2835),log(23590375)-log(167382319104), log(7.35)-9*log(10), log(5.64)-15*log(10), log(4.42)-23*log(10),log(2.6)-33*log(10), log(8.78)-46*log(10)];
d_list_birkhoff=2:min(d_max,length(birkhoff_exact_volumes)+1);
dim_power=-(d_list_birkhoff-1).^2.*log(d_list_birkhoff);
birkhoff_exact_volumes=birkhoff_exact_volumes(d_list_birkhoff-1)+dim_power;



vol=zeros(10,length(d_list));
for i=1:length(d_list)
   d=d_list(i);
   mu=ones(d,1)/d;
   nu=ones(d,1)/d;
   for j=1:reps
        vol(j,i)=volume_transport_polytopes(mu,nu,N, burn_in,subspace_normalization,return_samples, W,[]);  
   end
    
end



ind= vol==-inf;
vol(ind)=nan;
means_mh=mean(vol,'omitnan');
sd_mh=std(vol,'omitnan');

scatter(d_list_birkhoff,birkhoff_exact_volumes,20,'g')
hold on;
scatter(d_list,means_mh,40,'k');
hold on;
scatter(d_list,means_mh-2*sd_mh,20,'m');
hold on;
scatter(d_list,means_mh+2*sd_mh,20,'m');
hold on;

legend('exact volumes', 'mean marginal', 'standard deviation')
title('volume of the birkhoff polytope')
xlabel('dimension')
ylabel('log volume')