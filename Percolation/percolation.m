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

figure()
title('B&W image of percolation system Cutoff p=0.6')
imshow(z)

[lw,num] = bwlabel(z,4);  % lw -mx of labels for each cluster. Each cluster gets a number.
                          % num is the number of clusters

figure()
title('Color labeled clusters. Cutoff p=0.6')
img = label2rgb(lw,'jet','k','shuffle');  % create colour image
imshow(img);  % show.

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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% More plotting: finding beta.


figure()
subplot(2,1,1);
hold all
subplot(2,1,2);
hold all
% from now on we will use another pc:
pc = 0.59275;

betas = zeros(length(lstart:lend),1);
stdev = zeros(length(betas),1);
lnP = [];
lnp = [];
dPdp_mx = [];
for lcount = lstart:lend
    j = lcount - lstart + 1; % iteration 1,2,3....
    ln_P = [];
    ln_p_pc = [];
    pc_index = 0;
    for i=1:length(p)
        if p(i) > pc
            ln_p_pc(end + 1) = log(p(i)-pc);
            ln_P(end+1) = log(P(i,lcount));
        end
        if p(i) < (pc + (p(2)-p(1))) 
            pc_index = i;
        end
    end

    dPdp = diff(ln_P)./diff(ln_p_pc);
    beta = mean(dPdp);

    stdev(j) = std(dPdp,1);
    betas(j) = beta; 

%diff = length(ln_P) - length(ln_p_pc)
%should_be_41 = length(ln_P(diff:end))
%add = should_be_41 - diff

%len_P = length(ln_P)
%len_p = length(ln_p_pc)
%ln_P
%ln_p_pc
%len_dPdp = length(dPdp)
%len_p = length(p(pc_index:end-1))


    subplot(2,1,1)
    plot(ln_p_pc,ln_P,'-o')
    title('finding \beta ')
    ylabel('log(P(p,L))')
    xlabel('log(p-pc)')
    set(gca,'FontSize',fontsize)
    drawnow

    subplot(2,1,2)
    plot(p(pc_index:end-1),dPdp,'-d')
    set(gca,'FontSize',fontsize)
    ylabel('d(log(P))/dlog(p)')
    xlabel('p')
    drawnow

    lnP(:,end+1) = ln_P(:);
    lnp(:,end+1) = ln_p_pc(:);
    dPdp_mx(:,end+1) = dPdp(:);

end

subplot(2,1,1);
title('loglog plot')
axis([-10 -1 -2 0])
legend(legends,'Location','NorthWest')
subplot(2,1,2);
Title = ['Estimate on \beta *** p_c = ' num2str(pc,'%.5f')];
title(Title)

hold off


% plot the best loglog plot, which is the second last one according to
% output values in the standard deviation stdev.

% temporarly, it is the last loglog plot that is being replotted...
figure()
plot(p(pc_index:end-1),dPdp_mx(:,(lend - lstart -1)),'m-d')
set(gca,'FontSize',fontsize);
xlabel('p')
ylabel('d(ln(P))/d(ln(p))')
Title = ['Estimate on \beta *** \beta = ' num2str(betas(end-1),'%.2f') ' for L = ' num2str(2^(lend-1),'%g') ];
title(Title)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% see how well beta fits the P(p,L) distribution


figure()
subplot(1,1,1);
hold all
legends1 = {};
legends2 = {};
len_betas = length(betas);
%stdev = zeros(len_betas,1);
for i=1:len_betas
    lcount = lstart + i - 1;
    beta = betas(i);
    beta_func = (p(:)-pc).^(beta);
    
    plot(p(:),beta_func,'-d')
    %hold on
    plot(p(:),P(:,lcount),'-o')
    
    %set(gca,'FontSize',fontsize)
    ylabel('P(p,L)')
    xlabel('p')
    alegend1 = ['(p-pc)^{\beta} ; \beta = ' num2str(beta,'%.3f')];
    alegend2 = ['P(p,L=' num2str(2^lcount,'%g') ')'];
    legends1{end+1} = alegend1;
    legends1{end+1} = alegend2;
    drawnow
    
    %mean_dev = mean(P(:,lcount)-beta_func(:));
    %mean_sqrt = mean(sqrt(P(:,lcount)-beta_func(:)));
    
    %stdev(i) = sqrt(mean_sqrt/len_betas - (mean_dev/len_betas)^2);
end

subplot(1,1,1);
title('Similarity between P(p,L) and (p-pc)^{\beta}')
%C = [legends1,legends2];
legend(legends1,'Location','SouthEast')
hold off

% standard deviations:


stdev






end