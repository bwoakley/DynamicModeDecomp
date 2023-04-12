clear all;
close all;
clc;
%format long;


L = 50;               % Length of signal
T = 1.5;              % Sample Period

t = (0:L-1)*T/L;        % Time vector


% figure;
% plot(t)

S = 0.7*cos(2*pi*t) + sin(2*pi*2*t);

figure;
plot(S)

Y = fft(S);
% abs(Y)


if false
    rescaledY = Y/L;
    
    realCoeff = real(rescaledY);
    imagCoeff = imag(rescaledY);
    
    figure;
    plot(realCoeff)
    hold on;
    plot(imagCoeff)
    hold off;
    legend('real','imag')
    title('raw freq')
end

kk = 1;     %Use top kk frequencies (2*kk+1 total modes) to approx Eval1
if 2*kk+1 > L
    disp('kk too large')
end
shiftY = fftshift(Y);
% abs(shiftY)

if mod(L,2) == 0
    middle = ceil(L/2) + 1;
else 
    middle = ceil(L/2);
end



range = middle - kk: middle + kk ;

truncY = shiftY( range );

shiftTruncY = ifftshift(truncY);

ifft(shiftTruncY)






% k = (0:L-1)/1.5;
% 
% figure;
% plot(k, realCoeff)
% hold on;
% plot(k,imagCoeff)
% hold off;
% legend('real','imag')
% title('rescaled freq')

% rescaledY


if false
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


end

