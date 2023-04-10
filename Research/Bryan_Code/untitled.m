clear all;
close all;
clc;
format long;


Fs = 1000;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = 1500;             % Length of signal
t = (0:L-1)*T;        % Time vector

S = 0.7*cos(2*pi*50*t) + sin(2*pi*20*t);

figure;
plot(S)

Y = fft(S);
rescaledY = Y/L;

realCoeff = real(rescaledY);
imagCoeff = imag(rescaledY);





k = (0:L-1)/1.5;
% k = (0:L/2)/1.5;

figure;
plot(k, realCoeff)
hold on;
plot(k,imagCoeff)
hold off;
legend('real','imag')




% rescaledY









kk = 1;     %Use kk modes to approx Eval1
approxEval1 = 0;
for i = 1:kk
    approxEval1 = approxEval1 + 0;
end

% approxEval1 = ifft(Y, 1500);
approxEval1 = ifft(Y, 100,'symmetric');

figure;
plot( real(approxEval1) )
hold on;
plot( imag(approxEval1) )
legend('real','imag')


% figure;
% plot( abs(S - approxEval1) )
% title('error')




