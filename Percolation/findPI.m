

%
% Program to generate P(p,L) and Pi(p,L)

nsample = 10;                  % number of grid samples        
p = (0.35:0.01:1.0);           % probability span
nx = size(p,2);                % number of probabilities to run for
lstart = 3;                    % system size start
lend = 7;                      % system size stop

Pi = zeros(nx,lend);           % Pi(p,L) = p(probabiliity of having percolation)           
P = zeros(nx,lend);            % P(p,L)  = p(a sight is set as 1) = true if z < p(i) 
lvalue = zeros(lend);    
clf reset;        

legends = zeros(1,(lend-lstart));                    % list of legends

for lcount = lstart:lend
    count = lcount
    lx = 2^count;                                    % make lx*ly area that has size (2^count)^2
    ly = lx;
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
    
    % plotting:
    
    legends((lcount - lstart + 1)) = count; %num2str(count,'%g');
    title('Probability of percolation as a function of the cutoff probability')
    subplot(2,1,1);
    plt1 = plot(p,P(:,lcount));
    xlabel('p'); ylabel('P(p,L)');
    set(gca,'FontSize',18)
    hold on
  
    title('Probability of a sight being set as a function of the cutoff probability')
    subplot(2,1,2);
    plt2 = plot(p,Pi(:,lcount));
    xlabel('p'); ylabel('Pi(p,L)');
    set(gca,'FontSize',18)
    drawnow
end

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

