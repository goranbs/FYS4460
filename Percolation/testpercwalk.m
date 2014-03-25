%
% testpercwalk.m
%
% Generate spanning cluster (l-r spanning)
lx = 100;
ly = 100;
p = 0.59274;
nstep = 1e5;
nnstep = nstep + 1;
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
imagesc(zz),axis equal ,axis tight
rz = 1.0*zz;
n = 1;
while (n<=1)
r = rand(nnstep ,1);
[w,n] = percwalk(rz,r,0);
end
x = w(1,:);
y = w(2,:);
hold on,plot(y,x);
hold off