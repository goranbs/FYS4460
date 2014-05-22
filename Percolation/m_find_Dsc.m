% FYS4460 - Unsorted systems and percolation - Percolation project 3.
% exercise m)
%             Find the masss Msc fo teh singly connected bonds as a
%             function of the system size L for p = pc and use this to
%             estimate the exponent Dsc: Msc \propto L^Lsc. Also find the
%             behaviour of Psc = Msc/(L^d) as a function of (p-pc)
%
% We use the function exwalk.m to find the mass Msc of the singly connected
% bonds as a function of the system size L.


clear all 
close all

%L = [10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125];
L = [];
N = 100;   % 10 + 5*50 = 260;
for i = 1:N
    L(end+1) = 50 + 5*i;
end

len_L = length(L);
Msc = zeros(1,len_L);

for i = 1:len_L
    Msc(i) = exwalk(L(i),0);
    i
end

log_Msc = log(Msc);
log_L = log(L);

% slope Dsc:
[p,x] = polyfit(log_L,log_Msc,1);
f = polyval(p,log_L);
stdev =  mean((log_Msc - f).^2);

% plotting
fontsize = 18;

Title = ['Msc(L) Dsc = ' num2str(p(1), '%.3f') ' \pm ' num2str(stdev,'%.4f')];
name = 'm_Msc.png';
h = figure();
plot(log_L,log_Msc,'b-*')
hold on
plot(log_L,f,'r-')
set(gca,'Fontsize',fontsize)
title(Title)
xlabel('log(L)')
ylabel('log(Msc)')
print(h,'-dpng',name)

