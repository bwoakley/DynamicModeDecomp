clear all;
close all;
clc;
%format long;


L = 100;               % Length of signal
% T = 1.5;              % Sample Period

% t = (0:L-1)*T/L;        % Time vector
t = (1:L)/L;        % Let's rescale time to be just one time unit...


% figure;
% plot(t)

S = 0.7*cos(2*pi*.5*t) + sin(2*pi*4*t);

kk = 2;     %Use top kk frequencies to approx Eval1. Note kk should be < L/2


Y = fft(S);

rescaledY = Y/L;

% realCoeff = real(rescaledY);
% imagCoeff = imag(rescaledY);

% figure;
% plot(realCoeff)
% hold on;
% plot(imagCoeff)
% hold off;
% legend('real','imag')
% title('raw freq')

Yabs = abs(rescaledY); %Just looking at modulus. Which frequencies have most energy?

last = ceil(L/2); %Last positive freq
   
Ypos = Yabs(1:last);

[sort_Ypos, idx_Ypos] = sort(Ypos,'descend');   %idx_Ypos contains the indices's with most energy
top_freq = idx_Ypos - 1;




shiftY = fftshift(Y);

if mod(L,2) == 0
    middle = ceil(L/2) + 1;
else 
    middle = ceil(L/2);
end





if kk > last
    disp('kk too large')
end

Sapprox = zeros(1,L);

for j = 1:kk

    freq = top_freq(j)
    Sapprox = Sapprox +  shiftY(middle+freq)*exp((2*pi*1i)*freq*t) + shiftY(middle-freq)*exp((-2*pi*1i)*freq*t);
  

end

% shiftY
% 
% Sapprox

Sapprox = Sapprox/L;

tempSapprox = Sapprox;     %For some reason, Sapprox is shifted one to the left. Let's shift it back...
Sapprox(1) = S(1);
for j = 2:L
    Sapprox(j) = tempSapprox(j-1);
end
        
% Sapprox

figure; 
plot(t, S)
hold on;
plot(t, Sapprox,'--')
hold off;
legend('approx','orig')


if false
    kk = 20;     %Use top kk frequencies (2*kk+1 total modes) to approx Eval1
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
    
    Ytilde = ifft(shiftTruncY);
    
    figure;
    plot(Ytilde)

end


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

