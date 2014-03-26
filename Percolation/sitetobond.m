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
end