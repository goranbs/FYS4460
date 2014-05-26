%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%               PROJECT 4 - FYS4460
%
% Find the diffusion constant on a 2D grid


clear all
close all

L = 30;
Lx = L;
Ly = L;
time = 3000;   % timesteps
Nwalks = 2000;  % number of random walks on perc-cluster.

pc = 0.59275;

p = 0.4434;
%p = 0.467;
%p = 0.335156;
%p = pc;
%p=0.7;
cutoff_prob = p;

z = rand(Lx,Ly) < p;


% test if we have percolation:
[lw,num] = bwlabel(z,4);
[val,perc] = perc_test(lw,Lx,Ly);

counter = 0;
while (val == 0)
    z = rand(Lx,Ly) < p;
    [lw,num] = bwlabel(z,4);
    [val,perc] = perc_test(lw,Lx,Ly);
    counter = counter + 1;
end

num_attempts = counter;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% random walk on percolation cluster

exp_R = zeros(time,1);
exp_n = 0;
exp_s = 0;
exp_e = 0;
exp_w = 0;
exp_z = 0;
num_attempts_find_perc_cluster = zeros(Nwalks,1);

for nwalks = 1:Nwalks
    
%Walker_nr = nwalks

% 1) find starting point:

rw = [1+ (round((Lx-1)*rand(1))), 1+(round((Ly-1)*rand(1)))];

counter2 = 0;
truth = false;  % default

for i = 1:length(perc)
    if (perc(i) == lw(rw(1), rw(2)))
        truth = true;
    end
end


while (truth == false)
    
    for i = 1:length(perc)
        if (perc(i) == lw(rw(1), rw(2)))
            truth = true;
        else
            rw = [1+ (round((Lx-1)*rand(1))), 1+(round((Ly-1)*rand(1)))];
        end
    end
    counter2 = counter2 + 1;
end


num_attempts_find_perc_cluster(nwalks,1) = counter2;
random_walk_start = rw;
perc_cluster = lw(rw(1),rw(2));
perc_clusters = perc;

% 2) start walking:

r0 = rw; % initial position in cluster
r = r0;  % position of walker


R = zeros(time,1);
dirs = zeros(time,1); 

% time for the big walk on the percolation cluster!
for t = 2:time
   [r,dir] = walk(z,Lx,Ly,r);
   dirs(t-1,1) = dir;
   R(t) = ((r(1) - r0(1))^2 + (r(2) -r0(2))^2); % displacement!
   exp_R(t) = exp_R(t) + R(t);
end



north = (length(dirs(dirs==2)));
south = (length(dirs(dirs==1)));
east = (length(dirs(dirs==4)));
west = (length(dirs(dirs==3)));
stayed = length(dirs(dirs==0));

exp_n = exp_n + north;
exp_s = exp_s + south;
exp_e = exp_e + east;
exp_w = exp_w + west;
exp_z = exp_z + stayed;
 
Walker = ['Done with walker nr ' num2str(nwalks,'%g')]
end

exp_R = exp_R./Nwalks;
avg = 'average number of walks:'
exp_n = exp_n/Nwalks
exp_s = exp_s/Nwalks
exp_e = exp_e/Nwalks
exp_w = exp_w/Nwalks
exp_z = exp_z/Nwalks


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% done with the walking, now its time to calculate the diffusion

%per now, we only øook at one random walker, andt then:
T = linspace(0,time-1,time);
D = (diff(exp_R)./6);             % do not need to divide by ./diff(T), because diff(T) = [1 1 1.....1]
D2 = exp_R./(6*T');               % do not need to divide here either...

[p,S] = polyfit(T(500:end),exp_R(500:end)',1);
[Slope,st] = polyval(p,T(500:end),S);

DD = p(1)/6.0; % Slope

ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);
stdev1 = ste(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
fontsize = 18;

Diffusion_cons = D2(end);
stdev = std(D(500:end) - D2(501:end) );

Title1 = ['Diffusion. D=' num2str(Diffusion_cons,'%.3f') ' \pm ' num2str(stdev,'%.4f')];
Title2 = ['Displacement. p=' num2str(round(cutoff_prob*100),'%02g') '%. Nwalks=' num2str(Nwalks,'%g') ' timesteps=' num2str(time,'%g')];
Title3 = ['Diffusion. D=' num2str(DD,'%.4f') ' \pm ' num2str(stdev1,'%4f')];
name1 = ['diffusion_perc1_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name2 = ['displacement_perc_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name3 = ['diffusion_perc2_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];

h1 = figure();
set(gca,'Fontsize',fontsize)
plot(T(1:end-1),D','b-')
hold on
plot(T,D2','r-')
title(Title1)
xlabel('time [s]');ylabel('D')
legend('derivative','approx')
hold off
print(h1,'-dpng',name1)

h2 = figure();
set(gca,'Fontsize',fontsize)
plot(T,exp_R)
hold on
plot(T(500:end),Slope,'r-*')
title(Title2)
xlabel('time'); ylabel('<R^2>')
legend('<R^2>','ax+b')
print(h2,'-dpng',name2)

h3 = figure();
set(gca,'Fontsize',fontsize)
plot(T,D2','r-')
title(Title3)
xlabel('time [s]');ylabel('D')
legend('D')
print(h3,'-dpng',name3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some outputs
num_attempts_generate_perc_cluster = num_attempts
avg_num_attempts_find_perc_cluster = mean(num_attempts_find_perc_cluster)
Diffusion_constant = DD
standard_deviation = stdev1


