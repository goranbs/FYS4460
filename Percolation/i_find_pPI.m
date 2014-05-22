% Percolation project 3 - FYS4460
%
% i) Find pPI for PI = x for L=25,50,100,200,400,800 lattice.
%    plot pPI as a function of L.

L = [25,50,100,200,400,800];  % lattice sizes
x = [0.8,0.3];                % PI = x, probability of having percolation

Nsamples = 10000;               % number of samples to run for
pmin = 0.35;                   % lowest cutoff value
pmax = 0.75;                   % largest cutoff value
p = [pmin:0.005:pmax];         % cutoffs
len_p = length(p);          
len_L = length(L);
PI = zeros(len_p,len_L);      % holds PI as function of cutoffvalues for different lattice sizes


for Lsize = 1:len_L
    l = L(Lsize);
    for nsample = 1:Nsamples
        z=rand(l,l);                                 % (L*L) matrix of random distributed numbers in [0,1]
        for i = 1:len_p
            zz = z<p(i);                             % binary mx. 
            [lw,num]=bwlabel(zz ,4);                 % lw - mx of clusters. num = N_total
            perc_y = intersect(lw(:,1),lw(:,l));     % perc_y = #spanning clusters in y direction
            perc_x = intersect(lw(1,:),lw(l,:));     % perc_x = #spanning clusters in x direction
            perc_u = union(perc_x ,perc_y);          % perc = #spanning clusters
            perc = find(perc_u > 0);
     
            if (~isempty(perc))                      % If we have percolation:
                PI(i,Lsize) = PI(i,Lsize) + 1;       %  

            end
        end
    end
end

PI = PI./Nsamples;

% Test plotting to see if we get the right dataset 
% plot(p,PI(:,1))
% hold on
% plot(p,PI(:,2))
% plot(p,PI(:,3))
% plot(p,PI(:,4))
% plot(p,PI(:,5))
% plot(p,PI(:,6))
% hold off

% save the data set for plotting:

PI = [PI p(:)];

filename = 'PI_lattices.dat';
fileID = fopen(filename,'w');
fprintf(fileID,'%58s \n', '\\PI as function of cutoff probability. For latticesizes L');
fprintf(fileID,'%g %g %g %g %g %g\n',L(1),L(2),L(3),L(4),L(5),L(6));
dlmwrite(filename,PI,'-append', 'delimiter', ' ', 'precision', 13)

