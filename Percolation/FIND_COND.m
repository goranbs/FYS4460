
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

end
