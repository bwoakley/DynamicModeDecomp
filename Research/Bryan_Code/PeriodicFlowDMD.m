clear all; 
%clc; 
close all;
% Flow setups

Case=2; %1: Gyre; 2: Cylinder

% % Quad Gyre
if Case==1
dim1=101;dim2=dim1;
x=linspace(0,2*pi,dim1);y=x;
[Y X]=meshgrid(y,x);
om=1*pi;
k=0;
for t=0:.05:20
    k=k+1;
%    Psi=sin(X+t*om).*sin(Y);
%    Psi=sin(X+sin(t*om)).*sin(Y);
    Psi=sin(X).*sin(Y)*sin(om*t);
    PsiT(:,k)=Psi(:);
end

elseif Case==2
% Cylinder
dim1=251;dim2=101;
x=linspace(-2,8,dim1);
y=linspace(-2,2,dim2);
[Y X]=meshgrid(y,x);
a=1;alpha=2;R0=.35;L=2;y0=.3;U=14;W=24;
s=1-exp(-(X-1).^2/alpha^2-Y.^2);
k=0;
for t=0:.1:20
    k=k+1;
    x1=1+mod(t,1)*L;
    x2=1+mod(t-.5,1)*L;
    g1=exp(-R0*((X-x1).^2+alpha^2*(Y-y0).^2));
    g2=exp(-R0*((X-x2).^2+alpha^2*(Y+y0).^2));
    h1=abs(sin(pi*t));h2=abs(sin(pi*(t-.5)));
    g=-W*h1.*g1+W*h2.*g2+U*Y.*s;
    f=1-exp(-a*((X.^2+Y.^2).^1/2-1).^2);
    Psi=f.*g;
    ind=find(X.^2+Y.^2<1);Psi(ind)=0;
    PsiT(:,k)=Psi(:);
%    contour(X,Y,Psi,[-20:20],'k');daspect([1 1 1]);drawnow
end
else%1D test
    dim1=101;dim2=1;
    x=linspace(-1,1,dim1);
    f=x.^2-3*x;
    k=0;
    for t=0:.1:20
        k=k+1;
        Psi=f*sin(t);
        PsiT(:,k)=Psi(:);
    end
end

%%
n=40;
X = PsiT(:,1:n);
X2 = PsiT(:,2:n+1);
[U0,S0,V0] = svd(X,'econ');
ss=diag(S0);ind=find(ss>ss(1)*1e-10);
ss(ind)
%%  Compute DMD (Phi are eigenvectors)
r = length(ind);  % truncate at r modes, which is singular value >1e-10 max
U = U0(:,1:r);
S = S0(1:r,1:r);
V = V0(:,1:r);
Atilde = U'*X2*V*inv(S);
[W,eigs] = eig(Atilde);
Phi = X2*V*inv(S)*W;

%% Reconstructing DMD
Nf=121;%number of time steps forward from second to last time step
x1=X(:,end);%Take the last snap shot from X
b0=Phi\x1;%Obtain initial condition from last snapshot
temporal=zeros(r,Nf);
mu=diag(eigs);%eigenvalue
for i=1:Nf
    temporal(:,i)=b0.*mu.^i;%Raise by eigenvalue to power of timestep
end
X_dmd=Phi*temporal;

%% Plot reconstruction
close all;
tt=5;%time steps forward from n-1 to observe

theta=0:.1*pi:2*pi;xx=cos(theta);yy=sin(theta);
figure1=figure('Position',[150 200 1200 300]);
axes1 = axes('Parent',figure1,'Position',[.07 .1 .4 .8]);
xx1=reshape(X_dmd(:,tt),dim1,dim2);
xx2=reshape(PsiT(:,n+tt),dim1,dim2);
if Case ==3
    plot(x,real(xx1));
else
    pcolor(x,y,real(xx1)');shading interp;colorbar;daspect([1 1 1])%forecasting step
end
hold on;
fill(xx,yy,'k')
title('Forecasted field at T=4.5','fontsize',18)
% axes2 = axes('Parent',figure1,'Position',[.55 .55 .4 .4]);
% if Case ==3
%     plot(x,xx2)
% else
%     pcolor(x,y,xx2');shading interp;colorbar;daspect([1 1 1])%true data
% end
% hold on;
% fill(xx,yy,'k')
axes3 = axes('Parent',figure1,'Position',[.55 .1 .4 .8]);
if Case ==3
    plot(x,real(xx1)'-xx2');
else
    pcolor(x,y,real(xx1)'-xx2');shading interp;colorbar;daspect([1 1 1])%true data
end
hold on;
fill(xx,yy,'k')
title('Error at T=4.5','fontsize',18)

% axes4 = axes('Parent',figure1,'Position',[.55 .1 .4 .4]);
% xx3=reshape(PsiT(:,n+1),dim1,dim2);
% if Case ==3
%     plot(x,real(xx1)'-xx3');
% else
%     pcolor(x,y,real(xx1)'-xx3');shading interp;colorbar;daspect([1 1 1])%true data
% end
% hold on;
% fill(xx,yy,'k')



%% Inner Product
xx1V=xx1(:)-mean(xx1(:));
xx2V=xx2(:)-mean(xx2(:));
InnerProduct=dot(xx1V,xx2V);
CosineAngleBetween= dot(xx1V,xx2V)/(norm(xx1V,'fro')*norm(xx2V,'fro'));



%% Plot DMD modes
theta=0:.01*pi:2*pi;xx=cos(theta);yy=sin(theta);

figure2=figure('Position',[150 200 1400 300]);
subplot(1,3,1)
semilogy([1:10],ss(1:10))
set(gca,'position',[0.07 .2 .14 .6])
Ut=reshape(Phi(:,1),size(Psi,1),size(Psi,2));Ut=Ut';
title('Singular values','fontsize',18)
xlabel('Modal number');ylabel('Singular value')
set(gca,'fontsize',18)
subplot(1,3,2)
pcolor(x,y,Psi');shading interp;daspect([1 1 1]);colorbar;hold on;
fill(xx,yy,'k')

set(gca,'position',[0.26 .1 .3 .8])
title('Stream function','fontsize',18)
xlabel('X');ylabel('Y')
set(gca,'fontsize',18)

Ut=reshape(Phi(:,2),size(Psi,1),size(Psi,2));Ut=Ut';
subplot(1,3,3)
pcolor(x,y,real(Ut));shading interp;daspect([1 1 1]);colorbar;hold on;
fill(xx,yy,'k')
xlabel('X');ylabel('Y')
set(gca,'fontsize',18)

set(gca,'position',[0.64 .1 .3 .8])
title('DMD mode 2','fontsize',18)

% for i=1:2
%     Ut=reshape(Phi(:,i),size(Psi,1),size(Psi,2));Ut=Ut';
%     subplot(1,3,i+1)
%     pcolor(real(Ut));shading interp;daspect([1 1 1]);
% %     subplot(1,2,2)
% %     pcolor(imag(Ut));shading interp;daspect([1 1 1]);
% %     pause(1)
% end

%%  Plot DMD spectrum
figure
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);

%% Plot mu
figure;
plot(real(mu),'*-')
hold on; 
plot(imag(mu),'o-')
hold off;
title('Real and imaginary parts of the eigenvalues $\mu$ of $\widetilde{A}$','interpreter','latex')
legend('Real part','Imaginary part')


% figure;
% plot(ss)
