function [Phi, lambda, b, Xdmd, S, Atilde] = DMD(X1,X2,pred,r_orig)

% This code is borrowed from KutzBrunton2016 book.

% function [Phi ,omega ,lambda ,b,Xdmd ] = DMD(X1,X2,r,dt)
% Computes the Dynamic Mode Decomposition of X1, X2
%
% INPUTS:
% X1 = X, data matrix
% X2 = X’, shifted data matrix
% Columns of X1 and X2 are state snapshots
% r = target rank of SVD
% dt = time step advancing X1 to X2 (X to X’)
%
% OUTPUTS:
% Phi , the DMD modes
% omega , the continuous-time DMD eigenvalues
% lambda , the discrete -time DMD eigenvalues
% b, a vector of magnitudes of modes Phi
% Xdmd, the data matrix reconstructed by Phi , omega , b

% DMD
[W, S, V] = svd(X1, 'econ');

adjustV = true;
if adjustV
    sizeX1 = size(X1);
    for j = 1:sizeX1(2)    %adjusting vectors in V,W so that they are always sampled from the right half of R^M. This should stop the sign switching of Atilde
        if V(1,j) < 0
            V(:,j) = -1*V(:,j);
            W(:,j) = -1*W(:,j);
        end
    end
else
    disp('consider adjusting V')
end
%errorSignSwitch = sum(sum(abs(X1-W*S*V')))  %Notice that adjusting V and W still gives a valid SVD

% for j = 1:5
%     headW(j,:) = W(j,:);
% end
% headW


ss=diag(S);
ind=find(ss>ss(1)*1e-10);
% ss(ind)

%%  Compute DMD (Phi are eigenvectors)
max_r = length(ind);  % truncate at r modes, which is singular value >1e-10 max
r = min([r_orig, size(W,2), max_r]);
if r < r_orig
    disp(['Singular value(s) < ss(1)*1e-10 detected, truncating from r=', num2str(r_orig), ' to r=', num2str(r)])
end

W_r = W(:, 1:r); % truncate to rank-r
S_r = S(1:r, 1:r);
V_r = V(:, 1:r);
Atilde = W_r' * X2 * V_r / S_r; % low -rank dynamics

[Z_r , D] = eig(Atilde);
Phi = X2 * V_r / S_r * Z_r; % DMD modes
lambda = diag(D); % discrete -time eigenvalues

%omega = log(lambda)/dt; % continuous-time eigenvalues

% Compute DMD mode amplitudes b
%x1 = X1(:, end);
%b = Phi\x1;
% DMD reconstruction
%mm1 = size(X1, 2); % mm1 = m - 1
%time_dynamics = zeros(r, mm1);
%t = (0:mm1 -1)*dt; % time vector
%for iter = 1:mm1 ,
%time_dynamics (:,iter )=(b.*exp(omega*t(iter )));
%end;
%Xdmd = Phi * time_dynamics ;


% Reconstructing DMD
Nf=pred;%number of time steps forward 
x1=X2(:,end);%Take the last snap shot from X2
b=Phi\x1;%Obtain initial condition from last snapshot
temporal=zeros(r,Nf+1);
for i=1:Nf+1
    temporal(:,i)=b.*lambda.^(i-1);%Raise by eigenvalue to power of timestep
end
Xdmd=Phi*temporal;









