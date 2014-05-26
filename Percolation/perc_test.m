function [val,perc] = perc_test(z,lx,ly)
val = 0;
% Test if z has a spanning percolation cluster
% return:
% val = 1    : true, we have percolation
% val = 0    : false, we do not have percolation
%
% We also need to return which label(s) are the percolation cluster(s)!
% perc holds the labels of the percolation clusters! :-)

perc_x = intersect(z(1,:),z(lx,:));
perc_y = intersect(z(:,1),z(:,ly));
perc = union(perc_x,perc_y);
perc(perc==0) = [];                    % delete perc == 0

if(length(perc) > 0)
    val = 1;
end
end