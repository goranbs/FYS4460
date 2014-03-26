function g = coltomat(z,x,y)
% Convert z(x*y) into a matrix of z(x,y)
% Transform this onto a nx x ny lattice

g = zeros(x,y);
for iy = 1:y
    i = (iy -1)*x + 1;
    ii = i + x - 1;
    g(:,iy) = z(i:ii);
end
