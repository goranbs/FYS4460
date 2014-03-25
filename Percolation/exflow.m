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
zfn = fn>limit;
zbb = (zzz + 2*zfn);
zbb = zbb/max(max(zbb));
subplot(2,2,4), imagesc(zbb);
title('BB and DE');
axis equal
%
% Written by Marin Soreng
% (C) 2004
%
%Calculates the effective flow conductance Ceff of the
%lattice A as well as the pressure P in every site.
function [P, Ceff] = FIND_COND(A, X, Y)
P_in = 1;
P_out = 0;
%Calls MK_EQSYSTEM.
[B C] = MK_EQSYSTEM(A, X, Y);
%Kirchhoff 's equations solve for P
P = B\C;
%The pressure at the external sites is added
%(Boundary conditions)
P = [P_in*ones(X, 1); P; P_out*ones(X, 1)];
%Calculate Ceff
Ceff = (P(end -2*X+1:end-X)-P_out)'*A(end -2*X+1:end-X,2)/(P_in -P_out);

%
% Written by Marin S reng
% (C) 2004
%
% Sets up Kirchoff 's equations for the 2D lattice A.
% A has X*Y rows and 2 columns. The rows indicate the site ,
% the first column the bond perpendicular to the flow direction
% and the second column the bond parallel to the flow direction.
%
% The return values are [B, t] where B*x = C. This is solved
% for the site pressure by x = B\C.
function [B, C] = MK_EQSYSTEM(A, X, Y)
% Total no of internal lattice sites
sites = X*(Y-2);
%Allocate space for the nonzero upper diagonals
main_diag = zeros(sites , 1);
upper_diag1 = zeros(sites -1, 1);
upper_diag2 = zeros(sites -X, 1);
%Calculates the nonzero upper diagonals
main_diag = A(X+1:X*(Y-1), 1) + A(X+1:X*(Y-1), 2) + A(1:X*(Y-2), 2) ...
+ A(X:X*(Y-1)-1, 1);
upper_diag1 = A(X+1:X*(Y-1)-1, 1);
upper_diag2 = A(X+1:X*(Y-2), 2);
main_diag(find(main_diag==0)) = 1;
%Constructing B which is symmetric , lower=upper diagonals.
B = sparse(sites , sites);
% B*u = t
B = - spdiags(upper_diag1 ,-1, sites , sites);
B = B + - spdiags(upper_diag2 ,-X, sites , sites);
B = B + B' + spdiags(main_diag , 0, sites , sites);
%Constructing C
C = sparse(sites , 1);
C(1:X) = A(1:X, 2);
C(end-X+1:end) = 0*A((end -2*X+1:end-X), 2);
function g = sitetobond(z)
%
% Function to convert the site network z(L,L) into a (L*L,2) bond
% network
% g(i,1) gives bond perpendicular to direction of flow
% g(i,2) gives bond parallel to direction of flow
% z(nx,ny) -> g(nx*ny ,2)
%
nx = size(z,1);
ny = size(z,2);
N = nx*ny;
%g = zeros(N,2);
gg_r = zeros(nx,ny);
% First , find these
gg_d = zeros(nx,ny);
% First , find these
gg_r(:,1:ny -1) = z(:,1:ny -1).*z(:,2:ny);
gg_r(:,ny) = z(:,ny);
gg_d(1:nx-1,:) = z(1:nx -1,:).*z(2:nx ,:);
gg_d(nx ,:) = 0;
% Then , concatenate gg onto g
ii = 1:nx*ny;
g = zeros(nx*ny ,2);
g(:,1) = gg_d(ii)';
g(:,2) = gg_r(ii)';
function g = coltomat(z,x,y)
% Convert z(x*y) into a matrix of z(x,y)
% Transform this onto a nx x ny lattice
g = zeros(x,y);
for iy = 1:y
i = (iy -1)*x + 1;
ii = i + x - 1;
g(:,iy) = z(i:ii);
end
