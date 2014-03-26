

%
% Program to generate P(p,L) and Pi(p,L)
nsample = 5;
p = (0.35:0.01:1.0);
nx = size(p,2);
lstart = 5;
lend = 7;
Pi = zeros(nx,lend);
P = zeros(nx,lend);
lvalue = zeros(lend);
clf reset;
for lcount = lstart:lend
    lx = 2^lcount;
    ly = lx;
    ll = lx*ly;
    for ns = 1:nsample
        z=rand(lx,ly);
        for i = 1:nx
            zz = z<p(i);
            [lw,num]=bwlabel(zz ,4);
            perc_y = intersect(lw(:,1),lw(:,ly));
            perc_x = intersect(lw(1,:),lw(lx ,:));
            perc_u = union(perc_x ,perc_y);
            perc = find(perc_u >0);
            if (length(perc)>0)
                Pi(i,lcount) = Pi(i,lcount) + 1;
                s = regionprops(lw,'Area');
                area = cat(1,s.Area);
                ar = sum(area(perc_u(perc)));
                P(i,lcount) = P(i,lcount) + ar/ll;
            end
        end
    end
    
    P(:,lcount) = P(:,lcount)/nsample;
    Pi(:,lcount) = Pi(:,lcount)/nsample;
    subplot(2,1,1);
    plot(p,P(:,lcount));
    xlabel('p'); ylabel('P(p,L)');
    hold on
    subplot(2,1,2);
    plot(p,Pi(:,lcount));
    xlabel('p'); ylabel('Pi(p,L)');
    drawnow
end
hold off

