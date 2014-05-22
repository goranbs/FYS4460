%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%               PROJECT 4 - FYS4460
%
% Find the diffusion constant on a 2D grid

L = 10;
Lx = L;
Ly = L;

p = 0.3;

z = rand(Lx,Ly) < p;

% test if we have percolation:
[lw,num] = bwlabel(z,4);
val = perc_test(lw,Lx,Ly);

counter = 0;
while (val == 0)
    z = rand(Lx,Ly) < p;
    [lw,num] = bwlabel(z,4);
    val = perc_test(lw,Lx,Ly);
    counter = counter + 1;
end

counter

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% random walk on percolation cluster

