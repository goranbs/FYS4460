% finding an exponent.

close all
clear all

x = linspace(0,10,100);
x_squared = x.*x;


ln_x = log(x);
ln_x_squared = log(x_squared);

figure()
plot(ln_x,ln_x_squared)
xlabel('log(x)')
ylabel('log(x^2)')

ddx = diff(ln_x_squared)./diff(ln_x);


len_ddx =length(ddx)
len_x =length(x)

figure()
plot(x(1:end-1),ddx)
xlabel('x')
ylabel('ddx')


rectangularity = 1.0;             % cubic = 1
L = 100;                           % system size
r = rand(L,rectangularity*L);      % system
p = 0.6;                           % cutoff value. 60% of the values generated in r < p.
z = r <p;                          % binary matrix.

figure()
imshow(z)

[lw,num] = bwlabel(z,4);  % lw -mx of labels for each cluster. Each cluster gets a number.
                          % num is the number of clusters

figure()
img = label2rgb(lw,'jet','k','shuffle');  % create colour image
imshow(img);  % show.
