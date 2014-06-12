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
Nwalks = 300;  % number of random walks on perc-cluster.

pc = 0.59275;

%p = 0.335;    % part1
%p = 0.379;    % part2. run3
%p = 0.4434;   % part1.
%p = 0.45;      % part2. extra! 
%p = 0.467;    % part1.
%p = pc;
%p = 0.501;    %part2. run2
%p = 0.664;      % part2. run1, run4
p=0.7;
%p = 0.9;
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
    if (counter > 5000)
        counter = counter
        num_attempts = 'this does not work well!'
        break;
    end
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

%----------------------------- WALKERS -----------------------------------
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

%---------------------------- TIME-LOOP ----------------------------------
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
D = (diff(exp_R)./4);             % do not need to divide by ./diff(T), because diff(T) = [1 1 1.....1]
D2 = exp_R./(4*T');               % do not need to divide here either...

[p,S] = polyfit(T(500:end),exp_R(500:end)',1);
[Slope,st] = polyval(p,T(500:end),S);

DD = p(1)/4.0; % Slope

ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);
stdev1 = ste(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
fontsize = 18;

Diffusion_cons = D2(end);
stdev = std(D(500:end) - D2(501:end) );

%------------------------------- NAMES and TITLES ------------------------
Title1 = ['Diffusion. D=' num2str(Diffusion_cons,'%.3f') ' \pm ' num2str(stdev,'%.4f')];
Title2 = ['Displacement. p=' num2str(round(cutoff_prob*100),'%02g') '%. Nwalks=' num2str(Nwalks,'%g') ' timesteps=' num2str(time,'%g')];
Title3 = ['Diffusion. D=' num2str(DD,'%.4f') ' \pm ' num2str(stdev1,'%4f')];
Title4 = ['D=' num2str(DD,'%.4f') ' , p=' num2str(cutoff_prob,'%.2f') ' , Nwalks=' num2str(Nwalks,'%g')];
Title5 = ['D=' num2str(DD,'%.4f') ' , p=' num2str(cutoff_prob,'%.2f') ' , Nwalks=' num2str(Nwalks,'%g')];
name1 = ['diffusion_perc1_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name2 = ['displacement_perc_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name3 = ['diffusion_perc2_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name4 = ['log_log_diffusion_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];
name5 = ['log_log2_diffusion_p' num2str(round(cutoff_prob*100),'%02g') '_.png'];

%------------------------------PLOTTING-----------------------------------
h1 = figure();
set(gca,'Fontsize',fontsize)
plot(T(1:end-1),D','b-')
hold on
plot(T,D2','r-')
title(Title1)
xlabel('time [s]');ylabel('D')
leg1 = legend('derivative', 'approx');
set(leg1,'Location','SouthEast')
hold off
print(h1,'-dpng',name1)

h2 = figure();
set(gca,'Fontsize',fontsize)
plot(T,exp_R)
hold on
plot(T(500:end),Slope,'r-*')
title(Title2)
xlabel('time'); ylabel('<R^2>')
leg2 = legend('<R^2>', 'ax+b');
set(leg2,'Location','SouthEast')
print(h2,'-dpng',name2)

h3 = figure();
set(gca,'Fontsize',fontsize)
plot(T,D2','r-')
title(Title3)
xlabel('time [s]');ylabel('D')
legend('D')
print(h3,'-dpng',name3)

%--------------------------------------------------

ppc = (cutoff_prob - pc);
log_R = log(exp_R);
log_T = log(T);
log_Tppc = log_T + log(ppc);

END = round(time/30);

[p1,mu1] = polyfit((log_T(2:END))',log_R(2:END),1);
[p2,mu2] = polyfit((log_Tppc(5:END))',log_R(5:END),1);
f1 = polyval(p1,(log_T(2:END))');
f2 = polyval(p2,(log_Tppc(5:END))');

%mubeta = mu - beta;
%k_marked = nu/(2*nu + mu - beta);
%--------------------------------------------------

h4 = figure();
set(gca,'FontSize',fontsize)
plot(log_T,log_R,'r-')
hold on 
plot(log_T(2:END),f1,'b-*')
title(Title4)
leg4 = legend('log(<r^2>)','polyfit');
set(leg4,'Location','SouthEast')
xlabel('log(t)');ylabel('log(<r^2>)')
print(h4,'-dpng',name4)
hold off


h5 = figure();
set(gca,'FontSize',fontsize)
plot(log_Tppc,log_R,'r-')
hold on
plot(log_Tppc(5:END),f2,'b-*')
title(Title5)
leg5 = legend('log(<r^2>)','polyfit');
set(leg5,'Location','SouthEast')
xlabel('log(t(p-p_c))');ylabel('log(<r^2>)')
print(h5,'-dpng',name5)
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some outputs
num_attempts_generate_perc_cluster = num_attempts
avg_num_attempts_find_perc_cluster = mean(num_attempts_find_perc_cluster)
Diffusion_constant = DD
standard_deviation = stdev1

k = p1(1)/2
mubeta = p2(1)


