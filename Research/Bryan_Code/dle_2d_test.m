X=xnow;Y=ynow;
%X=analyticSolnX; Y=analyticSolnY; 

xss=xa;yss=ya;


%X(1,255)
%X(1,256)
%X(128,1)

%[G12 G11] = gradient(X,yss,xss);
[G11, G12] = gradient(X,xss,yss);


[G21, G22] = gradient(Y,xss,yss);

%figure;
%quiver(xa,ya,G11,G12)
%title('Green1')

%G11(1,256)
%G11(128,1)

%G12(1,1)
%G12(128,1)

%figure;
%quiver(xa,ya,G22,G21)
%title('quiver2')


%compute Cauchy-Green tensor
CG11 = G11.*G11 + G21.*G21 ;
CG12 = G11.*G12 + G21.*G22 ;

CG21 = G12.*G11 + G22.*G21 ;
CG22 = G12.*G12 + G22.*G22 ;

Tr=CG11+CG22;
Del=CG11.*CG22-CG12.*CG21;
lam=Tr/2+sqrt(Tr.^2-4*Del)/2;
dle=log(lam)/2;
[XSS, YSS]=meshgrid(xss,yss);

% scrsz = get(0,'ScreenSize');
% figure1=figure('Position',[50 50 scrsz(3)/3 scrsz(4)/1.5]);
% % subplot(3,1,1)
figure;
pcolor(XSS,YSS,dle);shading interp; colorbar; daspect([1 1 1]); 
title('FTLE')

% set(gca,'fontsize',24)
% daspect([1 1 1]); 

FTLE(:,:,ti)=dle;
% eval(ssss);

%figure;
%pcolor(XSS,YSS,X);shading interp; colorbar; daspect([1 1 1]); 
%title('X')

% figure;
% pcolor(XSS,YSS,G21);shading interp; colorbar; daspect([1 1 1]); 
% title('G11')
% 
% figure;
% pcolor(XSS,YSS,G22);shading interp; colorbar; daspect([1 1 1]); 
% title('G12')
% 
% min(min(G21))
% max(max(G21))