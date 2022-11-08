%X=xnow;Y=ynow; %Use these to compute net difference
xss=xa;yss=ya;
X = xafter; Y = yafter; %Use these to compute incremental diff
%xss = xbefore; yss = ybefore;

[G11, G12] = gradient(X,xss,yss);
[G21, G22] = gradient(Y,xss,yss);

%compute Cauchy-Green tensor
CG11 = G11.*G11 + G21.*G21 ;
CG12 = G11.*G12 + G21.*G22 ;

CG21 = G12.*G11 + G22.*G21 ;
CG22 = G12.*G12 + G22.*G22 ;

Tr=CG11+CG22;
Del=CG11.*CG22-CG12.*CG21;
lam=Tr/2+sqrt(Tr.^2-4*Del)/2;
%dle=log(lam)/(2);  
%dle=log(lam)/(2*INTTIME);            %Use this for net difference
dle=log(lam)/(2*deltat);    %Use this for incremental diff
[XSS, YSS]=meshgrid(xss,yss);

% scrsz = get(0,'ScreenSize');
% figure1=figure('Position',[50 50 scrsz(3)/3 scrsz(4)/1.5]);
% % subplot(3,1,1)


% figure;
% pcolor(dle);shading interp; colorbar; daspect([1 1 1]); 
% title('FTLE')

% set(gca,'fontsize',24)
% daspect([1 1 1]); 

%FTLE(:,:,ti)=dle;

% eval(ssss);









