%
% exflow.m
%
clear all; clf;
% First , find the backbone
% Generate spanning cluster (l-r spanning)

lx = 10;
ly = 10;
p = 0.5927;
ncount = 0;
perc = [];
while (size(perc ,1)==0)
    ncount = ncount + 1;
    if (ncount >1000)
        return
    end
    z=rand(lx,ly)<p;
    [lw,num]=bwlabel(z,4);
    perc_x = intersect(lw(1,:),lw(lx ,:));
    perc = find(perc_x >0)
end

s = regionprops(lw,'Area');
clusterareas = cat(1,s.Area);
maxarea = max(clusterareas);
i = find(clusterareas==maxarea);
zz = lw == i;
% zz now contains the spanning cluster

% Transpose
zzz = zz';

% Generate bond lattice from this
g = sitetobond(zzz);

% Generate conductivity matrix
[p c_eff] = FIND_COND(g,lx,ly);

% Transform this onto a nx x ny lattice
x = coltomat(full(p),lx,ly);
P = x.*zzz;
g1 = g(:,1);
g2 = g(:,2);
z1 = coltomat(g1,lx,ly);
z2 = coltomat(g2,lx,ly);

% Plotting
subplot(2,2,1), imagesc(zzz);
title('Spanning cluster')
axis equal
subplot(2,2,2), imagesc(P);
title('Pressure');
axis equal

f2 = zeros(lx,ly);
for iy = 1:ly-1
    f2(:,iy) = (P(:,iy) - P(:,iy+1)).*z2(:,iy);
end
    f1 = zeros(lx,ly);
for ix = 1:lx-1
    f1(ix ,:) = (P(ix ,:) - P(ix+1,:)).*z1(ix ,:);
end

% Find the sum of absolute fluxes into each site
fn = zeros(lx,ly);
fn = fn + abs(f1);
fn = fn + abs(f2);
fn(:,2:ly) = fn(:,2:ly) + abs(f2(:,1:ly -1));
fn(:,1) = fn(:,1) + abs((P(:,1) - 1.0).*(zzz(:,1)));
fn(2:lx ,:) = fn(2:lx ,:) + abs(f1(1:lx -1,:));

subplot(2,2,3), imagesc(fn);
title('Flux');
axis equal

limit = 0.01;
zfn = fn>limit;
zbb = (zzz + 2*zfn);
zbb = zbb/max(max(zbb));
subplot(2,2,4), imagesc(zbb);
title('BB and DE');
axis equal

