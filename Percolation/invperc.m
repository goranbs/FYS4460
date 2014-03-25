%
% invperc.m
%
% Example program for studying invasion percolation problems
% NOTE: This is not an optimal but an educational algorithm
%
L = 100;
% system size
p = (0.0:0.01:0.7);
perc = 0;
% flag to signal if other end is reached
nbetween = 1;
nstep = 0;
nend = numel(p);
nstop = 0;
z = rand(L,L);
% Random distribution of thresholds
pcluster = zeros(L,L);
while ((nstop==0)&&(nstep <nend))
nstep = nstep + 1;
p0 = p(nstep);
zz = z<p0;
[lw,num] = bwlabel(zz ,4);
leftside = lw(:,1);
i = find(leftside >0);
leftnonzero = leftside(i);
uniqueleftside = unique(leftnonzero);
cluster = ismember(lw,uniqueleftside);
pcluster = pcluster + cluster;
if (mod(nstep ,nbetween)==0)
imagesc(pcluster),axis equal , axis tight , colorbar ,drawnow
end
% Check if it has reached the right hand side
rightside = lw(:,L);
ir = find(rightside >0);
rightnonzero = rightside(ir);
span = intersect(leftnonzero ,rightnonzero);
if (numel(span)>0)
nstop = 1;
% spanning
end
end
p0
imagesc(pcluster),axis equal , axis tight , colorbar ,drawnow