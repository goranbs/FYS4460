%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   PROJECT 4 - FYS4460
%
%  Make a percolation cluster in 3D. Make sure that we have percolation.
%  Drop a randomwalker into the percolation cluster. Set P(r,0) = delta(r)
%  Measure the distance the random walker travels as a function of time.
%  That is, P(r,t).
%  Diffusion can be measured from (1/2d)(d<r^2>/dt) = D.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function perc_proj_4()

% system properties:
L = 10;                    % system size
Lx = L;
Ly = L;
Lz = L;
V = Lx*Ly*Lz;              % Volume

%p = 0.335156;              % probability of having percolation
p = 0.7;

z = zeros(Lx,Ly,Lz);                 % cube of zeros.
z = generate_system(z,Lx,Ly,Lz,p);   % system

%truth = perc_test(z,Lx,Ly,Lz);
%if (truth == 0)
%    z = generate_system(z,Lx,Ly,Lz,p);
%end



for kk = 1:Lz
    z(:,:,kk)
end

% making

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Random walk on the percolation cluster



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions

    function [sys] = generate_system(sys,lx,ly,lz,prob)
        
    for i = 1:lx
        for j = 1:ly
            for k = 1:lz
                random = rand(1);
                sys(i,j,k) = random < prob;
            
            end
        end
    end
    
    end

    function [val] = perc_test(sys,lx,ly,lz)
        % test if the system has a percolation cluster
        perc_x = intersect(sys(1,:,:),sys(lx,:,:));
        perc_y = intersect(sys(:,1,:),sys(:,ly,:));
        perc_z = intersect(sys(:,:,1),sys(:,:,lz));
        perc = union(perc_x,perc_y,perc_z);
        perc(perc == 0) = [];
        
        val = 0;
        if (len(perc) > 0)
            val = 1;
        end
    end
end