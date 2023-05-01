clear all;
close all;
clc;
format long;


kk = 256;

dt = .02;                   %time step size

Lx=20;
Ly=10;
rate = 1;  %Exp rate of decay

xa=linspace(-Lx,Lx,kk);
ya=linspace(-Ly,Ly,kk);

[XST, YST]=meshgrid(xa,ya);

ii = 0;

%Polar coordinates with center at xc,yc
xc = 0;
yc = 0;
Xt = XST - xc;
Yt = YST - yc;
R = sqrt(Xt.^2 + Yt.^2);

Z = sin(2*pi*R).*exp(-rate*R);

surf(XST,YST,Z,'EdgeColor','none')
% surf(Z)





