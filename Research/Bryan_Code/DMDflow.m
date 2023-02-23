clear all;
close all;
clc;
format long;


pred = 90;                   %pred = number of time steps forward to predict.
N = 10;                      %How many state snapshots to use (the X1 and X2 will have N-1 columns)
if N+pred>100
    disp('Need more data')
end


%Activate block to plot error for various r
%Here I used:
% pred = 90;                   
% N = 10; 
if false

    %r = 2;   %How many modes to use
    figure;
    rVec = [1,2,3,4,5,6,7,8,9,10];
    %rVec = [9,10];
    
    for r = rVec
        
        [error, lambda] = DMDpred(pred,N,r);
        plot(error)
        hold on;
        
    end
    
    hold off;
    legend('r=1','r=2','r=3','r=4','r=5','r=6','r=7','r=8','r=9','r=10')
    title('Error for different r truncations, given N = 10.')
    
    figure;
    plot(real(lambda),'*-')
    hold on; 
    plot(imag(lambda),'o-')
    hold off;

end



r = 5;   %How many modes to use
figure;
    
[error, lambda,S] = DMDpred(pred,N,r);
plot(error)

figure;
plot(real(lambda),'*-')
hold on; 
plot(imag(lambda),'o-')
hold off;

figure;
plot(diag(S))









