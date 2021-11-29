rng('shuffle')
%for finding the volume of general transport polytopes

%the number of nxm contingency tables with specified row and column sums 
%can be used to approximate the volume of transport polytopes. Examples
%2,3,4 show that this is surprisingly accurate!
%we can turn this idea on it's head: we can use our volume computation
%algorithms to approximate the number of contingency tables with given row
%and column sums (currently not implemented)

N=10000;
burn_in=1000;
reps=10;
log_scale=true;
subspace_normalization=false;
return_samples=false;
log_scale=true;
W=@(m)ones(length(m),1);
node_num=10;

%First example: 2x2 transport polytope
p=.7;
q=.1;
mu=[p;1-p];
nu=[q;1-q];
exact_volume=log(min(p,q)-max(0,p+q-1));

volume_transport_polytopes_rep(mu,nu,N, burn_in,subspace_normalization,return_samples, W,reps,exact_volume,log_scale,node_num);


    




%second example: from Beck &Pixton 2005 https://arxiv.org/pdf/math/0202267.pdf
subspace_normalization=false;

mu=[3046, 5173, 6116, 10928]';
nu=[182,778,3635,9558,11110]';
exact_volume=log(23196436596128897574829611531938753);
s=sum(mu);
mu=mu/s;
nu=nu/s;
exact_volume=exact_volume-(length(mu)-1)*(length(nu)-1)*log(s);

volume_transport_polytopes_rep(mu,nu,N, burn_in,subspace_normalization,return_samples, W,reps,exact_volume,log_scale,node_num);


%third example: from Beck &Pixton 2005 https://arxiv.org/pdf/math/0202267.pdf
mu=[338106, 574203,678876,1213008]';
nu=[20202,142746,410755,1007773,1222717]';
exact_volume=log(316052820930116909459822049052149787748004963058022997262397);
s=sum(mu);
mu=mu/s;
nu=nu/s;
exact_volume=exact_volume-(length(mu)-1)*(length(nu)-1)*log(s);

volume_transport_polytopes_rep(mu,nu,N, burn_in,subspace_normalization,return_samples, W,reps,exact_volume,log_scale,node_num);

%fourth example: from Beck &Pixton 2005 https://arxiv.org/pdf/math/0202267.pdf
mu=[30201, 59791, 70017, 41731, 58270]';
nu=[81016, 68993, 47000, 43001, 20000]';
exact_volume=log(24640538268151981086397018033422264050757251133401758112509495633028);
s=sum(mu);
mu=mu/s;
nu=nu/s;
exact_volume=exact_volume-(length(mu)-1)*(length(nu)-1)*log(s);
volume_transport_polytopes_rep(mu,nu,N, burn_in,subspace_normalization,return_samples, W,reps,exact_volume,log_scale,node_num);

function[ mn, sd]=volume_transport_polytopes_rep(mu,nu,N, burn_in,subspace_normalization,return_samples, W,reps,exact_volume,log_scale,node_num)

vol=zeros(10,1);
for j=1:reps
    vol(j)=volume_transport_polytopes(mu,nu,N, burn_in,subspace_normalization,return_samples, W,node_num); 
end
ind=vol==-inf;
vol(ind)=nan;
if ~log_scale
    vol=exp(vol);
    exact_volume=exp(exact_volume);
end
mn=mean(vol,'omitnan');
sd=std(vol,'omitnan');
disp(strcat('mean:',num2str(mn)));
disp(strcat('standard deviation:',num2str(sd)));
disp(strcat('exact volume:',num2str(exact_volume)));

end