% g) Plotting n(s,p)s^tau vs. F(s/sxi) to estimate sxi 
%
% Read output data from ClusterNumDensity.txt 
% produced from the MATLAB script ClusterNumberDensity.m
% From this text file we find n(s,p) which we use to estimate
% sigma in sxi = |p-pc|^(1/sigma).

fileID = fopen('ClusterNumDensity.dat');
Nsp = textscan(fileID,'%f64 %f64 %f64 %f64 %f64 %f64','delimiter',' ' ,'commentStyle', '\\');
fclose(fileID);

nsp1 = Nsp{1};
nsp2 = Nsp{2};
nsp3 = Nsp{3};
nsp4 = Nsp{4};
nsp5 = Nsp{5};
s = Nsp{6};

p = [nsp1(1) nsp2(1) nsp3(1) nsp4(1) nsp5(1)];  % p

nsp1 = nsp1(2:end);                 
nsp2 = nsp2(2:end);
nsp3 = nsp3(2:end);
nsp4 = nsp4(2:end);
nsp5 = nsp5(2:end);

nsp = [nsp1 nsp2 nsp3 nsp4 nsp5];              % n(s,p)
s = s(2:end);                                  % s

tau = 2.06;                                    % from tablevalues
pc = 0.59275;                                  % critical cutoff prob.

sigma = 36/91.0;                               % tabulated sigma value

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fontsize = 18;
legends = {};
h = figure();
hold all

for i=1:length(p)
    
    F = abs(p(i)-pc)^(1/sigma)*s(:);
    K = nsp(:,i).*(s(:).^(tau));
    
    alegend = ['p=' num2str(p(i),'%.3f')];
    legends{end+1} = alegend;
    plot(log10(F),log10(K),'-')
    set(gca,'FontSize',fontsize)
    xlabel('log(s|p-p_c|^{1/\sigma})');ylabel('log(s^{\tau}n(s,p))');
    drawnow
end

Title = ['Estimating s_{\xi}'];
title(Title)
legend(legends,'Location','SouthWest')
hold off

    