function percolation()
% FYS4460-Unsorted Systems and Percolation
%
% Project III - Percolation

L = 100;            % system size
r = rand(L,L);      % system
p = 0.6;            % cutoff value. 60% of the values generated in r < p.
z = r <p;           % binary matrix.
[lw,num] = bwlabel(z,4);  % lw -array of labels for each connected cluster
                          % num is the number of clusters
                          
img = label2rgb(lw,'jet','k','shuffle');  % create colour image
image(img);  % show.

s1 = regionprops(lw,'Area');        % s= struct array with filds: Area
area = cat(1,s1.Area);              % area is the are of the labels.

s2 = regionprops(lw,'BoundingBox'); % s= struct array with fields: BoundingBox 
bbox = cat(1,s2.BoundingBox);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a) make a feature to calculate P(p,L) for various p.

total_area = L*L;
N_total = num;
%N_spanning = ?
%imgsize=size(img)     % L x L x 3
%bboxsize=size(bbox)   % num x 4
%areasize=size(area)   % num x 1

%C = intersect(A,B);    % What is common for A and B
%D = union(A,B);        % the union of A and B.


end