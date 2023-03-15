clear all;
close all;
clc;
format long;

noStep = 6;
Atilde =  0.997989959900209;


a = zeros(1,noStep);
b = zeros(1,noStep);
b(1) = 0.998;

for i = 1:noStep       
    a(i) = (1-.02*i/10);
    b(i+1) = b(i)*Atilde;
end

a

b

figure;
plot(a)
hold on;
plot(b,'o');
hold off;

