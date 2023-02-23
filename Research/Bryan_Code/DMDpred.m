function [error,lambda, S] = DMDpred(pred,N,r)
%%This function will return the predictions for various choices of 
% N,pred, and r

% clear all;
% close all;
% clc;
% format long;

% %Test:
% row1 = linspace(1,10,10);
% row2 = linspace(2,11,10);
% data  = [ row1; row2];
% X1 = data(:,1:9);
% X2 = data(:,2:10);
% [Phi ,omega ,lambda ,b,Xdmd] = DMD(X1,X2,2,1);
% D = diag(lambda);
% 
% %pred = number of time steps forward to predict.
% pred = 2;
% Xpred = zeros(2,pred+1);
% Xpred(:,1) = X2(:, size(X2,2) );   %Instantiate: first column of Xpred will be the last column of X2
% for i = 1:pred
%     current = Xpred(:,i);
%     pseudo = Phi\current;
%     Xpred(:,i+1) = Phi * D * pseudo;
% end



no_of_steps = N + pred;     %Total number of velocity data to aquire (also get data past N to compare to our predictions).
dt = .02;                   %time step size

kk = 256;
rows = (256^2)*2;
stateVecs = zeros(rows,no_of_steps);

for i = 1:no_of_steps

        %data file name
        ss='Turb';

        ii = i + 1000;
        st=strcat('../Cases/',ss,'/bin0',num2str(ii));
        fid=fopen(st,'rb');
        data=fread(fid,[1 1],'*float');
        data=fread(fid,[1 inf],'*double');
        fclose(fid);
        data=reshape(data,kk,kk,3);
        u1 = data(:,:,1)' ;
        v1 = data(:,:,2)' ;

%         figure;
%         quiver(u1,v1)

        tempU = reshape(u1,[],1);
        tempV = reshape(v1,[],1);
        temp = [tempU; tempV];

        stateVecs(:,i) = temp;

end


X1 = stateVecs(:,1:N-1);
X2 = stateVecs(:,2:N);

[Phi ,omega ,lambda ,b,Xdmd, S, uMax] = DMD(X1,X2,r,dt);
D = diag(lambda);


Xpred = zeros(rows,pred+1);
Xpred(:,1) = X2(:, size(X2,2) );   %Instantiate: first column of Xpred will be the last column of X2
for i = 1:pred
    current = Xpred(:,i);
    pseudo = Phi\current;
    Xpred(:,i+1) = Phi * D * pseudo;
end


    %Now analize the error at each time step
error = zeros(1,pred);
%errorC = zeros(1,pred);

for i = 1:pred+1

    temp = Xpred(:,i) - stateVecs(:,N-1+i) ;
    error(i) = sumabs(real(temp));
    %errorC(i) = sumabs(imag(temp));

end

error = error/(rows*uMax);     %Normalize it
error = error(2:end);   %The first entry is the last entry of X2, so let's drop it.
% figure;
% plot(error,'o-')
% hold on;
% plot(errorC,'*-')
% hold off;






