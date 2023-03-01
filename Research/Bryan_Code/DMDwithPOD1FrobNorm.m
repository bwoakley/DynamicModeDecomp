clear all; close all; %clc
load ../output/091226_DLEseriesB.mat

dim1=size(DLE,1);dim2=size(DLE,2);
Nt=120;
u=DLE(1:end,1:end,1:Nt+1);
for i=1:Nt+1
    tu=u(:,:,i);
    M(:,i)=tu(:);%-mean(tu(:));
end 

% M=M';
size(M);

for Ti=1:100
n=20;
X = M(:,Ti:n-1+Ti);
X2 = M(:,Ti+1:n+Ti);
[U,S,V] = svd(X,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 20;  % truncate at r modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde = U'*X2*V*inv(S);

XX = M(:,Ti+1:n+Ti);
XX2 = M(:,Ti+2:n+1+Ti);
[U,S,V] = svd(XX,'econ');

%%  Compute DMD (Phi are eigenvectors)
r = 20;  % truncate at r modes
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
Atilde1 = U'*XX2*V*inv(S);
errs(Ti)=norm(Atilde-Atilde1);
end

plot(1:100,errs)
xlabel('Time');title('$||\tilde{A}_{i+1}-\tilde{A}_i||$','Interpreter', 'latex');ylabel('error');
set(gca,'fontsize',18)