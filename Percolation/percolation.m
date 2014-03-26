function percolation()
% FYS4460-Unsorted Systems and Percolation
%
% Project III - Percolation

clf reset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) make a feature to calculate P(p,L) for various p.


%N_spanning = ?
%imgsize=size(img)     % L x L x 3
%bboxsize=size(bbox)   % num x 4
%areasize=size(area)   % num x 1

%C = intersect(A,B);    % What is common for A and B
%D = union(A,B);        % the union of A and B.

rectangularity = 1.0;             % cubic = 1
L = 100;                           % system size
r = rand(L,rectangularity*L);      % system
p = 0.6;                           % cutoff value. 60% of the values generated in r < p.
z = r <p;                          % binary matrix.
[lw,num] = bwlabel(z,4);  % lw -mx of labels for each cluster. Each cluster gets a number.
                          % num is the number of clusters

figure()
title('cutoff p=0.6')
img = label2rgb(lw,'jet','k','shuffle');  % create colour image
image(img);  % show.

s1 = regionprops(lw,'Area');        % s= struct array with filds: Area
area = cat(1,s1.Area);              % area is the are of the labels.

s2 = regionprops(lw,'BoundingBox'); % s= struct array with fields: BoundingBox 
bbox = cat(1,s2.BoundingBox);

total_area = L*L;
N_total = num;

nsample = 10;                  % number of grid samples        
p = (0.35:0.01:1.0);           % probability span
nx = size(p,2);                % number of probabilities to run for
lstart = 3;                    % system size start
lend = 7;                      % system size stop

Pi = zeros(nx,lend);           % Pi(p,L) = p(probabiliity of having percolation)           
P = zeros(nx,lend);            % P(p,L)  = p(a sight is set as 1) = true if z < p(i) 
lvalue = zeros(lend);    

%clf reset;        

legends = zeros(1,(lend-lstart));                    % list of legends
pc = 0;
counter = 0;
for lcount = lstart:lend
    count = lcount
    lx = 2^count;                                    % make lx*ly area that has size (2^count)^2
    ly = rectangularity*lx;
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
            if P(i,lcount) < 0.8
                pc = pc + p(i);
                counter = counter + 1;
            end
        end
    end
    
    % plotting:
    %Title = ['pc= ' num2str(pc,'%g')]
    legends((lcount - lstart + 1)) = count; %num2str(count,'%g');
    %title(Title)
    title('Probability of percolation as a function of the cutoff probability, and P(p,L)')
    subplot(2,1,1);
    plt1 = plot(p,P(:,lcount),'b-o');
    xlabel('p'); ylabel('P(p,L)');
    set(gca,'FontSize',18)
    hold on
  
    
    %title('Probability of a sight being set as a function of the cutoff probability')
    subplot(2,1,2);
    %title(Title)
    plt2 = plot(p,Pi(:,lcount),'b-o');
    xlabel('p'); ylabel('Pi(p,L)');
    set(gca,'FontSize',18)
    drawnow
end

pc = pc/counter
%legends
%legends(:) = num2str(legends(:),'%s')
%llength = length(legends);
%legneds2 = []
%for i = 1:llength
%    legends2(i) = legends(i);
%end

%legend(legends)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We see that the probability of a sight being set as a function of the
% cutoff probability is closing up on the value 0.6.
%
% 
hold off




end