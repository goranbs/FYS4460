function percolation()
% FYS4460-Unsorted Systems and Percolation
%
% Project III - Percolation


%N_spanning = ?
%imgsize=size(img)     % L x L x 3
%bboxsize=size(bbox)   % num x 4
%areasize=size(area)   % num x 1

%C = intersect(A,B);    % What is common for A and B
%D = union(A,B);        % the union of A and B.

clear all
close all
clf reset

fontsize = 18;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) make a feature to calculate P(p,L) for various p.


rectangularity = 1.0;             % cubic = 1
L = 100;                           % system size
r = rand(L,rectangularity*L);      % system
p = 0.6;                           % cutoff value. 60% of the values generated in r < p.
z = r <p;                          % binary matrix.
[lw,num] = bwlabel(z,4);  % lw -mx of labels for each cluster. Each cluster gets a number.
                          % num is the number of clusters

%figure()
title('cutoff p=0.6')
img = label2rgb(lw,'jet','k','shuffle');  % create colour image
%image(img);  % show.

s1 = regionprops(lw,'Area');        % s= struct array with filds: Area
area = cat(1,s1.Area);              % area is the are of the labels.

s2 = regionprops(lw,'BoundingBox'); % s= struct array with fields: BoundingBox 
bbox = cat(1,s2.BoundingBox);

total_area = L*L;
N_total = num;
p_c = 0.59275;                 % experimental value of the critical probability cutoff
p_min = p_c + 0.001;           % p_min > p_c
nsample = 10;                  % number of grid samples        
p = (p_min:0.01:1.0);           % probability span
nx = size(p,2);                % number of probabilities to run for
lstart = 3;                    % system size start
lend = 7;                      % system size stop

Pi = zeros(nx,lend);           % Pi(p,L) = p(probabiliity of having percolation)           
P = zeros(nx,lend);            % P(p,L)  = p(a sight is set as 1) = true if z < p(i) 
if p_min > p_c
    beta = zeros(nx,lend);         % exponential: P(p,L) ~ (p-pc)^beta.
end
lvalue = zeros(lend);    

%clf reset;        

%get(gca,'ColorOrder')
figure()
subplot(2,1,1);
hold all
subplot(2,1,2);
hold all

%legends = zeros(1,(lend-lstart));                    % list of legends
legends = {};
pc = 0;
counter = 0;
pc_list = [];

for lcount = lstart:lend
    lx = 2^lcount;                                   % make lx*ly area that has size (2^count)^2
    ly = rectangularity*lx;                          % I mean: lx*(rectangularity*ly).
    Total_area = lx*ly;
    for ns = 1:nsample
        z=rand(lx,ly);                               % (lx*ly) matrix of random distributed numbers in [0,1]
        for i = 1:nx
            zz = z<p(i);                             % binary mx. 
            [lw,num]=bwlabel(zz ,4);                 % lw - mx of clusters. num = N_total
            perc_y = intersect(lw(:,1),lw(:,ly));    % perc_y = #spanning clusters in y direction
            perc_x = intersect(lw(1,:),lw(lx ,:));   % perc_x = #spanning clusters in x direction
            perc_u = union(perc_x ,perc_y);          % perc = #spanning clusters
            perc = find(perc_u > 0);
            if (length(perc) > 0)                    % If we have percolation:
                Pi(i,lcount) = Pi(i,lcount) + 1;     %  
                s = regionprops(lw,'Area');          %
                area = cat(1,s.Area);                % 
                ar = sum(area(perc_u(perc)));        % find area, ar, of percolation 
                P(i,lcount) = P(i,lcount) + ar/Total_area;  % P(p,L)
            end
        end
    end
    
    P(:,lcount) = P(:,lcount)/nsample;     % mean
    Pi(:,lcount) = Pi(:,lcount)/nsample;   % mean
    beta(:,lcount) = log(P(:,lcount) - (p(:)-p_c));
    
    % find pc:
    for i = 1:nx
        if Pi(i,lcount) > 0.2
            if Pi(i,lcount) < 0.8
                pc = pc + p(i);
                counter = counter + 1;
                pc_list(counter) = p(i);
            end
        end
    end
    
    % plotting:
    %Title = ['pc= ' num2str(pc,'%g')]
    
    alegend = num2str(lx,'%g');
    alegend = ['L=' alegend];
    legends{end+1} = alegend;
    %legends(1,ncounts) = alegend
    
    subplot(2,1,1);
    plot(p,P(:,lcount),'-o');
    xlabel('p'); ylabel('P(p,L)');
    set(gca,'FontSize',fontsize)
    drawnow
  
    
    %title('Probability of a sight being set as a function of the cutoff probability')
    subplot(2,1,2);
    %title(Title)
    plot(p,Pi(:,lcount),'-o');
    xlabel('p'); ylabel('Pi(p,L)');
    set(gca,'FontSize',fontsize)
    drawnow
end

pc = pc/counter

subplot(2,1,1);
title('Probability of a site being within a spanning cluster as function of the cutoff probability p')
subplot(2,1,2);
Title = ['Probability of percolation. pc ~ ' num2str(pc,'%.2f')];
title(Title)
legend(legends,'Location', 'SouthEast')
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% More plotting: finding beta.

p2 = p(2:end) - p_c;
figure()
subplot(2,1,1);
hold all
subplot(2,1,2);
hold all
for lcount = lstart:lend

pfunc = (p2(:)).^beta(2:end,lcount);

%alegend2 = {'P(p,L) ' num2str(lcount) ,'(p-pc)^{\beta} ' num2str(lcount)};
%legends2(lcount - lstart +1) = alegend2;
subplot(2,1,1);
set(gca,'FontSize',fontsize)
plot(p(2:end),P(2:end,lcount),'-o')
xlabel('p'); ylabel('P(p,L)');
drawnow

subplot(2,1,2);
set(gca,'FontSize',fontsize)
plot(p(2:end),pfunc(:),'-*')
xlabel('p'); ylabel('(p-p_c)^{\beta}')
drawnow
end

subplot(2,1,1)
Titl = ['Probability of cite being within spaning cluster.'];
title(Titl)
%legend(legends2)

hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We see that the probability of a sight being set as a function of the
% cutoff probability is closing up on the value 0.6.
%



end