% Example of use of the walk routine
% Generate spanning cluster (l-r spanning)

function exwalk()
lx =64;
ly = 64;
p = 0.585;
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
imagesc(zz);

% Display spanning cluster
% Run walk on this cluster
[l,r] = walk(zz);
zzz = l.*r;

% Find points where both l and r are non-zero
zadd = zz + zzz;
subplot(2,2,1), imagesc(zz);
subplot(2,2,2), imagesc(zadd);
subplot(2,2,3), imagesc(zzz >0);
subplot(2,2,4), imagesc(l+r>0);
function [left ,right] = walk(z)
%
% Left turning walker
%
% Returns left: nr of times walker passes a site
%
% First , ensure that array only has one contact point at left and
% right end: topmost points chosen
%

nx = size(z,1);
ny = size(z,2);
i = find(z(1,:)>0);
iy0 = i(1);

% starting point for walker
ix0 = 1;

% stopping point for walker
% First do left -turning walker
dirs = zeros(4,2);
dirs(1,1) = -1;
dirs(1,2) = 0;
dirs(2,1) = 0;
dirs(2,2) = -1;
dirs(3,1) = 1;
dirs(3,2) = 0;
dirs(4,1) = 0;
dirs(4,2) = 1;
nwalk = 1;
ix = ix0;
iy = iy0;
dir = 1;

% 1=left , 2 = down , 3 = right , 4 = up;
left = zeros(nx,ny);
while (nwalk >0)
    left(ix,iy) = left(ix,iy) + 1;
    
    % Turn left until you find an occupied site
    nfound = 0;
    while (nfound==0)
        dir = dir - 1;
        if (dir <1)
            dir = dir + 4;
        end
        % Check this direction
        iix = ix + dirs(dir ,1);
        iiy = iy + dirs(dir ,2);
        if (iix==nx+1)
            nwalk = 0;
            % Walker escaped
        iix = nx;
        ix1 = ix;
        iy1 = iy;
        end
        
        % Is there a site here?
        if (iix >0)
            if (iiy >0)
                if (iiy<ny+1)
                    if (z(iix,iiy)>0)
                        % there is a site here , move here
                        ix = iix;
                        iy = iiy;
                        nfound = 1;
                        dir = dir + 2;
                        if (dir >4)
                            dir = dir - 4;
                        end
                    end
                end
            end
        end
    end
end

%left;
nwalk = 1;
ix = ix0;
iy = iy0;
dir = 1;
% 1=left , 2 = down , 3 = right , 4 = up;
right = zeros(nx,ny);
while (nwalk >0)
    right(ix,iy) = right(ix,iy) + 1;
    % ix,iy
    % Turn right until you find an occupied site
    nfound = 0;
    while (nfound==0)
        dir = dir + 1;
        if (dir >4)
            dir = dir - 4;
        end
        % Check this direction
        iix = ix + dirs(dir ,1);
        iiy = iy + dirs(dir ,2);
        if (iix==nx+1)
            if (iy==iy1)
                nwalk = 0;
                % Walker escaped
                iix = nx;
            end
        end
        % Is there a site here?
        if (iix >0)
            if (iiy >0)
                if (iiy<ny+1)
                    if (iix<nx+1)
                        if (z(iix,iiy)>0)
                            % there is a site here , move here
                            ix = iix;
                            iy = iiy;
                            nfound = 1;
                            dir = dir - 2;
                            if (dir <1)
                                dir = dir + 4;
                            end
                        end
                    end
                end
            end
        end
    end
end
end
end