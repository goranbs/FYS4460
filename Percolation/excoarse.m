
function excoarse()
%
% excoarse.m
%
% Example of use of the coarsening procedure
%
z = rand(512,512)<0.58;
%
% Set up array for f

f(1) = 0;
f(2) = 0;
f(3) = 0;
f(4) = 1;
f(5) = 0;
f(6) = 0;
f(7) = 0;
f(8) = 1;
f(9) = 0;
f(10) = 0;
f(11) = 0;
f(12) = 1;
f(13) = 1;
f(14) = 1;
f(15) = 1;
f(16) = 1;
[lz,nz] = bwlabel(z,4);
imgz = label2rgb(lz);
zz = coarse(z,f);
[lzz,nzz] = bwlabel(zz ,4);
imgzz = label2rgb(lzz);
zzz = coarse(zz,f);
[lzzz ,nzzz] = bwlabel(zzz ,4);
imgzzz = label2rgb(lzzz);
zzzz = coarse(zzz,f);
[lzzzz ,nzzzz] = bwlabel(zzzz ,4);
imgzzzz = label2rgb(lzzzz);
subplot(2,2,1), image(imgz);
axis equal
subplot(2,2,2), image(imgzz);
axis equal
subplot(2,2,3), image(imgzzz);
axis equal
subplot(2,2,4), image(imgzzzz);
axis equal

function zz = coarse(z,f)
    % The original array is z
    % The transfer function is f given as a vector with 16 possible places
    % f applied to a two-by-two matrix should return
    % the renormalized values
    %
    % The various values of f correspond to the following
    % configurations of the two-by-two region that is renormalized ,
    % where I have used X to mark a present site , and 0 to mark an
    % empty sites
    %
    % 1 00 5 00 9 00 13 00
    % 00 X0 0X XX
    %
    % 2 X0 6 X0 10 X0 14 X0
    % 00 X0 0X XX
    %
    % 3 0X 7 0X 11 0X 15 0X
    % 00 X0 0X XX
    %
    % 4 XX 8 XX 12 XX 16 XX
    % 00 X0 0X XX
    %
    nx = size(z,1);
    ny = size(z,2);
    if (mod(nx ,2)==1)
        return
    end
    if (mod(ny ,2)==1)
        return
    end
    nx2 = floor(nx/2);
    ny2 = floor(ny/2);
    zz = zeros(nx2,ny2);
    x=zeros(2,2);
    for iy = 1:2:ny
        for ix = 1:2:nx
            x = 1 + z(ix,iy)*1 + z(ix,iy+1)*2 + z(ix+1,iy)*4 + z(ix+1,iy+1)*8;
            xx = f(x);
            zz((ix+1)/2,(iy+1)/2) = xx;
        end
    end
end

end