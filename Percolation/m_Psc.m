% FYS4460 - Unsorted systems and percolation - Percolation project 3.
% exercise m) part 2
%
%            find the behaviour of Psc = Msc/(L^d) as a function of (p-pc)
%
% We use the function exwalk.m to find the mass Msc of the singly connected
% bonds as a function of (p-pc). We have then changed the functionality of
% exwalk so that it takes both L and p as variables.


clear all
close all

pc = 0.59275;
d = 2;  % 2 dimentions
L = [2*64 3*64 4*64];
pmin = 0.56;
pmax = 0.76;

%p = pmin:0.01:pmax;
%p = [0.56 0.57 0.58 0.59 pc 0.6 0.61 0.62 0.63 0.64];
p = [0.58 0.59 pc 0.6 0.61 0.62 0.63 0.64 0.66 0.7];

Nsamples = 10;

Psc = zeros(length(L),length(p));

for systemsize = 1:length(L)
    SysSize = L(systemsize)
    for prob = 1:length(p)
        Prob = p(prob)
        Msc = 0;
        for sample = 1:Nsamples
            Msc = Msc + exwalk(L(systemsize),p(prob),0);
        end
        Psc(systemsize,prob) = (double(Msc)/Nsamples)/(L(systemsize)^d);
    end
end

ppc = p-pc;

% plotting

fontsize = 18;
legends = {};
Title = ['Evolution of Psc as funciton of (p-p_c)'];
name = ['m_2_Psc.png'];

h = figure();
subplot(1,1,1);
hold all
for i = 1:length(L)
    plot(ppc, Psc(i,:))
    set(gca,'FontSize',fontsize)
    xlabel('(p-p_c)'),ylabel('Psc');
    legends{end+1} = ['L=' num2str(L(i),'%g')];
    drawnow
end

subplot(1,1,1)
title(Title)
legend(legends)
print(h,'-dpng',name)