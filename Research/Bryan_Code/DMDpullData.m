function [X1, X2, stateVecs] = DMDpullData( pred, N, flowCase, kk, start_index )
%%This function will return the data X1, X2 for various choices of 
% N, pred, and start_index

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

no_of_steps = N + pred;     %Total number of (velocity) data columns to aquire (also get data past N to compare to our predictions).
dt = .02;                   %time step size

rows = (kk^2)*2;
stateVecs = zeros(rows,no_of_steps);

for i = 1:no_of_steps

        if flowCase == 1   %then choose 'turb' flow
            
            %data file name
            ss='Turb';
    
            ii = i + 1000 + start_index-1+20;
            st=strcat('../Cases/',ss,'/bin0',num2str(ii));
            fid=fopen(st,'rb');
            data=fread(fid,[1 1],'*float');
            data=fread(fid,[1 inf],'*double');
            fclose(fid);
            data=reshape(data,kk,kk,3);
            u1 = data(:,:,1)' ;
            v1 = data(:,:,2)' ;
            omega = data(:,:,3)' ;

%             figure;
            pcolor(omega); shading interp; drawnow;

        elseif flowCase == 2

            %Length of domain
            L = 3;
            x=linspace(0,L,kk*3+1); 
            y=linspace(0,L,kk*3+1);
            xa=x(kk+1:2*kk);ya=y(kk+1:2*kk);
            [XST, YST]=meshgrid(xa,ya);
            u1 = -1*(XST-1.5);
            v1 = YST-1.5; 
        elseif flowCase == 3
  
            ii = i + start_index-1;
 
            %Length of domain
            L = 3;
            x=linspace(0,L,kk*3+1); 
            y=linspace(0,L,kk*3+1);
            xa=x(kk+1:2*kk);ya=y(kk+1:2*kk);
            [XST, YST]=meshgrid(xa,ya);
            u1 = -1*(1-dt*i/10)*(XST-1.5);
            v1 = (1-dt*i/10)*(YST-1.5);  

        end
%         figure;
%         quiver(u1,v1)

        tempU = reshape(u1,[],1);
        tempV = reshape(v1,[],1);
        temp = [tempU; tempV];

        stateVecs(:,i) = temp;

end


X1 = stateVecs(:,1:N-1);
X2 = stateVecs(:,2:N);

% D = diag(lambda);


% Xpred = zeros(rows,pred+1);
% Xpred(:,1) = X2(:, size(X2,2) );   %Instantiate: first column of Xpred will be the last column of X2
% for i = 1:pred
%     current = Xpred(:,i);
%     pseudo = Phi\current;
%     Xpred(:,i+1) = Phi * D * pseudo;
% end






