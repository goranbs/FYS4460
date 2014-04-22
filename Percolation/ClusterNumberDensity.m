function ClusterNumberDensity()
clear all
close all
fontsize = 18;

pc = 0.59275;             % for a given cutoff p;
L = 512;                  % system size 
rectangularity = 1.0;     % rectangularity
lx = L;                  
ly = rectangularity*lx;
Ld = lx*ly;               % system size
z = rand(lx,ly) < pc;

% find the distribution of cluster sizes in this system.
%
%       n(s,p) = N(s)/Ld

Nsamples = 50;
a = 1.2;
binmax = 10e5;
[lw,num] = bwlabel(z,4);
Cs = regionprops(lw,'Area');
areas = cat(1,Cs.Area);
%[xbin,dx,xbins] = logbin(areas,a,binmax);

% from below
p_min = 0.55;
p_max = 0.7;
%probs = p_min:0.01:(pc+0.1);  % from below pc -> pc
%probs = (pc+0.001):0.01:p_max;
probs = [0.55 0.57 pc 0.599 0.6]

h = figure();
subplot(1,1,1);
hold all
%subplot(2,1,2);
%hold all

legends = {};

n_sp     = [];   % cluster number density
alphas   = [];   % exponent
a_stdevs = [];   % standard deviations
a_means  = [];   % mean alpha value

N_bins = [];

for i = 1:length(probs)
    
    istart=i
    for sample = 1:Nsamples
        
        z = rand(lx,ly) < probs(i);
        
        [lw,num] = bwlabel(z,4);
        Cs = regionprops(lw,'Area');
        areas = cat(1,Cs.Area);
        
        % check if we have a percolating cluster.
        % if we do, then remove this from the custersizes Cs.
        
        perc_y = intersect(lw(1,:),lw(end,:));
        perc_x = intersect(lw(:,1),lw(:,end));
        perc_u = union(perc_y,perc_x);          % labels of percolating clusters
        perc_u(perc_u == 0) = [];               % perc_u == 0 returns indexes of perc_u = 0. Delete these indexes.
        areas(perc_u) = [];                     % Delete indexes in areas that is percolating clusters.
        
        [bins,bin_width,Nbins] = logbin(areas,a,binmax);
        
        if sample == 1
            nbins = Nbins(:);
        else
            nbins(:) = nbins(:) + Nbins(:);
        end
    end
    
    bins(:) = bins(:);
    nbins(:) = (nbins(:)./bin_width')/(Ld*Nsamples);
  
    N_bins(end+1,:) = nbins(:);
    if i==length(probs)
        N_bins(end+1,:) = bins(:);
    end

    ln_s = log10(bins);
    ln_n = log10(nbins);
%     
%     subplot(2,1,1)
%     plot(bins,nbins);
%     set(gca,'FontSize',fontsize)
%     xlabel('s');ylabel('n(s,p)')
%     drawnow


    subplot(1,1,1)
    if probs(i) == pc
        plot(ln_s,ln_n,'-')
        alegend = ['pc=' num2str(probs(i),'%.3f')];
        legends{end+1} = alegend;
    else
        plot(ln_s,ln_n,'--')
        
        alegend = ['p=' num2str(probs(i),'%.3f')];
        legends{end+1} = alegend;
    end
    set(gca,'FontSize',fontsize)
    xlabel('log(s)');ylabel('log(n(s,p))');
    drawnow
    
    % find the power-law exponent:
    
%     alpha = diff(ln_n)./(diff(ln_s));
%     
%     n_sp(end+1,:) = Nbins(:);
%     alphas(end+1,:) = alpha(:);
%     a_stdevs(end+1) = std(alpha);
%     a_means(end+1) = mean(alpha);
%    
    idone=i
end

%subplot(2,1,1);
%title('Cluster Number Density')

subplot(1,1,1);
Title = ['LogLog-plot *** Cluster Number Density. #samples=' num2str(Nsamples,'%g')];
title(Title)
legend(legends,'Location','NorthEast')
axis([0 6 -13 -2])
print(h,'-dpng','clusterNRdensity_loglog.png')

filename = 'ClusterNumDensity.dat';
fileID = fopen(filename,'w');
fprintf(fileID,'%87s \n', '\\Fist line probabilities: p. For column below p: n(s,p). Last cloumn: cluster sizes: s');
fprintf(fileID,'%.3f %.3f %.3f %.3f %.3f %g\n',probs(1),probs(2),probs(3),probs(4),probs(5),0);
dlmwrite(filename,N_bins','-append', 'delimiter', ' ', 'precision', 13)

% figure()
% hold all
% legends2 = {};
% %len_nsp = length(n_sp)
% %len_alphas = length(alphas)
% for i = 1:length(probs)
%    plot(n_sp(1:end,i),alphas(:,i),'-d')
%    xlabel('Cluster number density [A_{cluster}/L^2]');ylabel('\alpha');
%    
%    alegend = ['p' num2str(i,'%g') ' , \alpha =' num2str(a_means(i),'%.2f') ' \pm' num2str(a_stdevs(i),'%.3f')];
%    legends2{end+1} = alegend;
% end
% 
% title('Power-law exponent');
% legend(legends2)

end
