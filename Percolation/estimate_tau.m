
% n(s,pc) = C*s^(-tau)
clear all
close all
fontsize = 18;

pc = 0.59275;      % Critical cut-off probability.
Nsamples = 100;   % number of samples to get the statistics.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_min = 4;         % system sizes L = 2^k
k_max = 9;
a = 1.2;
binmax = 10e5;

legends = {};
h = figure(1);
hold all

log10_nsp = [];

for k = k_min:k_max

    L = 2^k
    for sample = 1:Nsamples
        
        z = rand(L,L) < pc;
        
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
        
        clear z;
    end
    
    bins(:) = bins(:);
    nbins(:) = (nbins(:)./bin_width')/(L*L*Nsamples);
    
    ln_s = log10(bins);
    ln_n = log10(nbins);
    
    ln_s(~isfinite(ln_n)) = [];
    ln_n(~isfinite(ln_n)) = [];
    len_n = length(ln_n);
    
    plot(ln_s(1:len_n),ln_n,'--')
    
    alegend = ['L=' num2str(L,'%g')];
    legends{end+1} = alegend;
    
    set(gca,'FontSize',fontsize)
    xlabel('log(s)');ylabel('log(n(s,p))');
    drawnow
    
    %log10_nsp(end+1,:) = ln_n(:);
    % dette fungerte dårlig....:
    if k == k_max
        ln_s(~isfinite(ln_n)) = [];
        ln_n(~isfinite(ln_n)) = [];
        len_n = length(ln_n);
        p = polyfit(ln_s(1:(end-4)), ln_n(1:(end-4))',1);
        f = polyval(p,ln_s);
        
        tau = -p(1);
        stdev = mean((ln_n' - f).^2);
        plot(ln_s,f,'-')
        set(gca,'FontSize',fontsize)
        xlabel('log(s)');ylabel('log(n(s,p))');
        drawnow
        alegend = ['\tau = ' num2str(tau,'%.3f') ' \pm ' num2str(stdev,'%.3f')];
        legends{end+1} = alegend;
        
        % or we could find the average slope for the 512x512 case...
%         
%         dnds = diff(ln_n)./diff(ln_s(1:len_n)');
%         dnds(~isfinite(dnds)) = [];             % delete inf-values.
%         
%         tau = -mean(dnds(1:end-2))
%         stdev = std(dnds(1:end-2),1)
%         
%         C = 0.05;
%         log_nsp = log10(C*bins.^(-tau));
%         
%         plot(log10(bins),log_nsp,'-')
%         xlabel('log(s)');ylabel('log(n(s,p))')
%         alegend = ['\tau = ' num2str(tau,'%.2f') ' \pm ' num2str(stdev,'%.3f')];
%         legends{end+1} = alegend;
%         drawnow
        
        log10_nsp(end+1,:) = ln_s;
        log10_nsp(end+1,:) = ln_n;
    end
    
    
    
end

Title = ['LogLog-plot *** Cluster Number Density. #samples= ' num2str(Nsamples,'%g')];
title(Title)
legend(legends,'Location','SouthWest')
print(h,'-dpng','Estimate_tau.png')                     % save plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write to file

filename = 'estimate_tau.txt';
fileID = fopen(filename,'w');
fprintf(fileID,'%33s %.3f %4f \n','\\ tau and its standard deviation',tau,stdev);
fprintf(fileID, '%12s %8s \n', '\\ log10_nsp','log10_s');
dlmwrite(filename, log10_nsp' , '-append', 'delimiter',' ', 'precision', 13)



