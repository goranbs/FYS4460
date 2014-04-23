% h) Mass scaling of percolation cluster.
%
% Find the mass M(L) for the percolating cluster

% fileID = fopen('estimate_tau.dat');
% C = textscan(fileID,'%f64 %f64', 'delimiter', ' ', 'commentStyle', '\\');
% fclose(fileID);
% 
% ln_n = C{1};  % ln_n for L=512 system
% ln_s = C{2};  % ln_s for L=512 system
% Ok, we can surely do this, but this is in fact not what we are asked!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We should plot M(L) for L = 2^k for k=4,5,....,11.



clear all
close all

fontsize = 18;
pc = 0.59275;      % Critical cut-off probability.
Nsamples = 200;    % number of samples to get the statistics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_min = 4;         % system sizes L = 2^k
k_max = 11;
a = 1.2;

Mass = [];
Syst = [];
for k = k_min:k_max

    L = 2^k
    M = 0;
    for sample = 1:Nsamples
        
        z = rand(L,L) < pc;
        
        [lw,num] = bwlabel(z,4);

        
        % check if we have a percolating cluster.
        % if we do, then let this be part of the statistics.
        % no other areas need to be stored.
        
        perc_y = intersect(lw(1,:),lw(end,:));
        perc_x = intersect(lw(:,1),lw(:,end));
        perc_u = union(perc_y,perc_x);          % labels of percolating clusters
        perc = find(perc_u > 0);                % indexes in perc_u larger than zero
        
        if (length(perc) > 0)                   % there exist other values than 0, we have percolation.
            Cs = regionprops(lw,'Area');
            areas = cat(1,Cs.Area);
            M = M + sum(areas(perc_u(perc)));   % add up the area/mass of the percolating clusers
        end
    end
    Mass(end+1) = M/Nsamples;
    Syst(end+1) = L;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the exponent

ln_M = log10(Mass);
ln_L = log10(Syst);

% microscopic way:
dMdL = diff(ln_M)./diff(ln_L);
D1 = mean(dMdL);
stdev1 = std(dMdL,1);

% macroscopic way:
p = polyfit(ln_L,ln_M,1);
f = polyval(p,ln_L);
D2 = p(1);
stdev2 = mean((ln_M -f).^2);

Title1 = ['Mass scaling of percolating cluster. D=' num2str(D1,'%.4f') ' \pm ' num2str(stdev1,'%.5f')];
filename1 = ['MasScaling1.png'];
Title2 = ['Mass scaling of percolating cluster. D=' num2str(D2,'%.4f') ' \pm ' num2str(stdev2,'%.5f')];
filename2 = ['MasScaling2.png'];

% Plotting

h = figure(1);
plot(ln_L,ln_M,'b-*')
set(gca,'FontSize',fontsize)
ylabel('log(M(L))');xlabel('log(L)');
legend('M(L)','Location','SouthEast');
title(Title1)
print(h,'-dpng',filename1)

h = figure(2);
plot(ln_L,ln_M,'b-*')
hold on
plot(ln_L,f,'g--')
set(gca,'FontSize',fontsize)
ylabel('log(M(L))');xlabel('log(L)');
legend('M(L)','L^D','Location','SouthEast')
title(Title2)
print(h,'-dpng',filename2)
hold off

D1
stdev1
D2
stdev2