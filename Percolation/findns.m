%
% findns.m
%
lx=256;
ly=lx;
ll=lx*ly;
nsample = 1;
logbinsize = 2;
logbinmax = ll;
p = (0.2:0.05:0.6);
nx = size(p,2);
Pi = zeros(nx ,1);
P = zeros(nx ,1);
for i = 1:nx
    for isample = 1:nsample
        z=rand(lx,ly);
        zz = z<p(i);
        [lw,num]=bwlabel(zz ,4);
        perc_y = intersect(lw(:,1),lw(:,ly));
        perc_x = intersect(lw(1,:),lw(lx ,:));
        perc_xy = union(perc_x ,perc_y);
        perc = find(perc_xy >0);
        s = regionprops(lw,'Area');
        clusterareas = cat(1,s.Area);
        if (length(perc)>0)
            % Set Pi
            Pi(i) = Pi(i) + 1;
            % Find P(p,L)
            ar = sum(clusterareas(perc_xy(perc)));
            P(i) = P(i) + ar/ll;
        end
        % Find the cluster number density , get rid of percolating clusters
        ind = (1:num);
        indnoP = setxor(ind,perc_xy(perc));
        % Do statistics on area(indnoP)
        clusta = clusterareas(indnoP);
        [x,dx,n] = logbin(clusta ,logbinsize ,logbinmax);
        if (isample==1)
            nnsp = n/ll;
            nnsp = nnsp'./dx;
            nsp = nnsp;
        else
            nnsp = n/ll;
            nnsp = nnsp'./dx;
            nsp = nsp + nnsp;
        end
    end

    P(i) = P(i)/nsample;
    Pi(i) = Pi(i)/nsample;
    nsp = nsp/nsample;
    ind2 = find(nsp >0);
    plot(log10(x(ind2)),log10(nsp(ind2)),'-o');
    xlabel('s');ylabel('n(s,p)');
    hold on; drawnow;
end
hold off;