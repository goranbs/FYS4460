function powelawdist()
% Determine the exponent of a power-law distribution.
clear all
close all

fontsize = 18;


% N = 1e4;
% nsamples = 10;
% p_min = 1.0e3;     % [0,1] -> [0,0.1] -> [0,0.01]
% p_max = 4.0e4;
% 
% 
% %z = rand(N,1).^(-3+1)
% 
% % first finding the cumulative distribution function P(Z>z).
% 
% Z = linspace(p_min,p_max,100); % values to check for.
% 
% len_z= length(Z);
% 
% P = zeros(len_z,1);     % Probability of z being less than Z.
% 
% for i = 1:len_z
%     
%     for k = 1:nsamples
%         
%         z = rand(N,1).^(-3+1);  % generate random distribution for z(i)
% 
%         for j = 1:length(z)
%             
%             if Z(i) > z(j)
%                 
%                 P(i,1) = P(i,1) + 1;
%             end
%         end
%    
%     end
%    
%     
% end
% 
% P(:,1) = P(:,1)/nsamples/N; 
% 
% figure()
% plot(Z,P)  % cumulative distribution
% 
% 
% %length(P)
% %length(z)
% %len_z
% 
% fZ = diff(P')./diff(Z); % dP/dz
% 
% %len_Z
% %len_fZ = length(fZ)
% 
% figure()
% plot(Z(1:end-1),fZ)     % c) find the distribution function
% 
% ln_Z = log(Z(:));
% ln_P = log(P(:));
% 

%-------------------------------------------------------------

n = 1e6;
kappa = 1-3;
z = rand(n,1).^(kappa);
zz = sort(z);
Nz = (1:n)/n;

h = figure(1);
plot(Nz,zz,'b-*')
set(gca,'FontSize',fontsize)
axis([0.999 1.0 0 1e10])
title('Distribtion of dataset. z= rand(n,1).^{(-3+1)}')
xlabel('Nz [1/n]'); ylabel('Sorted random value; z')
print(h,'-dpng','Z_power_law.png')

Nsamples = 10;

len_Z = 20;
Z_min = 0.0001;
Z_max = 1.0;
Z = linspace(Z_min,Z_max,len_Z);

P_Z = zeros(len_Z,1);

% for j = 1:len_Z
%     for sample = 1:Nsamples
%         z = rand(n,1).^(-kappa);
%         %zz = sort(z);
%         for i = 1:n
%             if Z > z(i)
%                 P_Z(j) = P_Z(j) + 1;
%             end
%         end
%     end
% end
% 
% 
% P_Z = P_Z(:)/(Nsamples*n);
% 
% figure()
% plot(Z,P_Z)

z = rand(n,1).^(kappa);
[bin,bin_width,n_values] = logbin(z,1.5,10^5);

Bin_width = bin_width;
Bins = bin;
N_values = n_values;

for samples = 1:(Nsamples-1)
    
    z = rand(n,1).^(kappa);
    [bin,bin_width,n_values] = logbin(z,1.5,10^5);
    
    Bins(:) = Bins(:) + bin(:);
    N_values(:) = N_values + N_values(:);
end

Bins(:) = Bins(:)/Nsamples;
N_values(:) = N_values(:)/Nsamples;

fZ = N_values./(Bin_width');
h = figure(2);
bar(Bins,fZ)
set(gca,'FontSize',fontsize)
axis([0 12e4 0 10])
xlabel('Size of number'); ylabel('Number of values present in dataset')
title('Bar plot')
print(h,'-dpng','Bar_plot_find_exponent.png')

% looks like a straight line to me :-)

% find the exponent:

alpha = diff(log(fZ))./(diff(log(bin))');
std_alpha = std(alpha,1)
alpha = mean(alpha)


h = figure(3);
plot(log(bin),log(fZ),'b-*')
set(gca,'FontSize',fontsize)
Title = ['LogLog-plot. \alpha = ' num2str(alpha,'%.3f') ' \pm ' num2str(std_alpha,'%.4f')];
title(Title)
xlabel('log(z)'); ylabel('log(f_Z)') 
print(h,'-dpng','LogLog_alpha_find_exponent.png')

h = figure(4);
plot(Bins, fZ)
set(gca,'FontSize',fontsize)
title('Probability density distribution')
xlabel('z'); ylabel('f_Z') 

print(h,'-dpng','prob_dens_find_exponent.png')

hold off

end