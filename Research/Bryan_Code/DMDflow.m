clear all;
close all;
clc;
format long;


pred = 90;                   %pred = number of time steps forward to predict.
N = 10;                      %How many state snapshots to use (the X1 and X2 will have N-1 columns)
if N+pred>100
    disp('Need more data')
end


%%Activate block to plot error for various r
%Here I used:
% pred = 90;                   
% N = 10; 
if false
    flowCase = 1; %flowCase decides what flow to use. flowCase = 1 means 'turb'
    figure;
    rVec = [1,2,3,4,5,6,7,8,9]; %Max r to take is N-1
    
    for r = rVec
        
        [error, lambda, S] = DMDpred(pred,N,r,flowCase);
        plot(error)
        hold on;
        
    end
    
    hold off;
    legend('r=1','r=2','r=3','r=4','r=5','r=6','r=7','r=8','r=9')
    title('Error for different r truncations, given N = 10.')
    xlabel('Time steps after N')
    ylabel('Error')


    figure;
    plot(real(lambda),'*-')
    hold on; 
    plot(imag(lambda),'o-')
    hold off;

end


%%Activate block to plot singular values and eigenvalues
%Here I used:
%N = 10
%r = 5, 6, or 9
if false

    flowCase = 1; %flowCase decides what flow to use. flowCase = 1 means 'turb'
    %N = 10;
    %pred = 100-N;
    r=9;
    %r = N-1;   %How many modes to use

%     figure;
%     [error, lambda,S] = DMDpred(pred,N,r,flowCase);
%     plot(error)
%     title('error')
    
    figure;
    plot(real(lambda),'*-')
    hold on; 
    plot(imag(lambda),'o-')
    hold off;
    title('Real and imaginary parts of the eigenvalues $\Lambda$ of $\widetilde{A}$','interpreter','latex')
    legend('Real part','Imaginary part')
    
    figure;
    plot(diag(S))
    title('Singular values of $X_1$','interpreter','latex')

end


%%Now fix r and vary N.
if true

    flowCase = 1; %flowCase decides what flow to use. flowCase = 1 means 'turb'
    r = 5;   %How many modes to use
    
    figure;
    nVec = [10,20,30,40,50]; %Max pred to take is 100-N
    
    for N = nVec
    
        pred = 100-N;
    
        [error, lambda, S] = DMDpred(pred,N,r,flowCase);
        plot(error)
        hold on;
        
    end
    
    hold off;
    legend('N=10','N=20','N=30','N=40','N=50')
    title('Error for different N, given r = 5.')
    xlabel('Number of time steps predicted')
    ylabel('Error')
end


%%Now try linear flow.
% ss='Line';
% kk = 256;
% 
% i=1;
% 
% ii = i + 1000;
% st=strcat('../Cases/',ss,'/bin0',num2str(ii));
% fid=fopen(st,'rb');
% data=fread(fid,[1 1],'*float');
% data=fread(fid,[1 inf],'*double');
% fclose(fid);
% data=reshape(data,kk,kk,3);
% u1 = data(:,:,1)' ;
% v1 = data(:,:,2)' ;
% 
% figure;
% quiver(u1,v1)


%Length of domain
L = 3;
x=linspace(0,L,256*3+1); 
y=linspace(0,L,256*3+1);
xa=x(257:512);ya=y(257:512);
[XST, YST]=meshgrid(xa,ya);
u1 = -1*(XST-1.5);
v1 = YST-1.5;
% figure;
% quiver(u1,v1)
