function [x,dx,n] = logbin(y,a,binmax)
%
% Returns a logarithmically binned dataset
% y is a vector of the data -set
% a is the bin size , that is, bins are from a^k to a^k+1
% x gives the centers of the bins
% dx gives the width of the bins
% n gives the number of points in each bin
%
% First , general a list of edges , smallest value of y is 1
%ymax = max(y);
ymax = binmax;
yedge = 1.0;
istep = 1;
yedgelast = 0;
while (yedgelast <=ymax)
    edge(istep) = yedge;    
    yyedge = floor(yedge*a);
    dy = yyedge - yedge;
    if (dy <=1.0)
        yyedge = yedge + 1.0;
    end
    yedgelast = yedge;
    yedge = yyedge;
    istep = istep + 1;
end

n = histc(y,edge);
dx = diff(edge);
nx = size(edge ,2);
x = 0.5*(edge(1:nx -1) + edge(2:nx));
n = n(1:nx -1);
end