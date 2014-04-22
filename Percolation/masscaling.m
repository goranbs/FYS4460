% h) Mass scaling of percolation cluster.
%
% Find the mass M(L) for the percolating cluster

fileID = fopen('estimate_tau.dat');
C = textscan(fileID,'%f64 %f64', 'delimiter', ' ', 'commentStyle', '\\');
fclose(fileID);

ln_n = C{1};  % ln_n for L=512 system
ln_s = C{2};  % ln_s for L=512 system
% Ok, we can surely do this, but this is in fact not what we are asked!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% We should plot M(L) for L = 2^k for k=4,5,....,11.






