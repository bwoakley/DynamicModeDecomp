close all
clear all

N = 20;
t = linspace(0,3*pi/2,N);

for i = 1:N-1

    v(i) = sin(t(i)) ;
    w(i) = sin(t(i+1)) ;
end

figure
plot(v,w)

v 

w