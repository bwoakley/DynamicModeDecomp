%% I want to plot the FTLE fields on the left and the moms on the right

% r_dir=' /media/philw/Thanos/fullcode_data/dg_data0/'; %Read the FTLE field files;
% rm_dir=' /media/philw/Thanos/fullcode_data/dg_data_test/'; %Where to read the moms from;
%     load_mom=strcat(' load',rm_dir,'G/full_moments.mat');
%     eval(load_mom);
    
    %xc=2*(50/128);
    %yc=2*(44/128);
    
    
    
     XLogan=linspace(0,4,256+1); XLogan=XLogan(1:end-1); XLogan=XLogan(:);
     YLogan=linspace(0,4,256+1); YLogan=YLogan(1:end-1); YLogan=YLogan(:);
%     
    [YLogan,XLogan]=meshgrid(YLogan,XLogan);
    XLogan=XLogan(:)-2;
    YLogan=YLogan(:)-2;
dt=.02;
kx=(.001^2)/2;


exponent = (.5*(XLogan).^2 + .5*(YLogan).^2);
%amplitude = 1 / (2 * sqrt(2*pi));  
% The above is very much different than Alan's "1./2*pi*sigma^2"
% which is the same as pi*Sigma^2 / 2.
%val       = amplitude  * 

     Ptot= 1 / (2*pi*10000)*exp(-exponent);
     RealP=Ptot/sum(Ptot(:));
     %break
for i=1:200
    
    ti=mod(i,20); ti(ti==0)=1;
    load_field=strcat(' load',r_dir,'ftle_field',num2str(ti),'.mat');
    eval(load_field);
    

    
    
    %Try reshaping the strain fields to find something...?
    
    
    
%      fieldx=cgStrainD(:,1).*cgStrainV(:,1)+cgStrainD(:,2).*cgStrainV(:,3);
%      fieldy=cgStrainD(:,1).*cgStrainV(:,2)+cgStrainD(:,2).*cgStrainV(:,4);
%      alpha=fieldx./(XLogan.*dt);
%      beta=fieldy./(YLogan.*dt);
%      Exx=-kx./alpha+(kx./alpha).*exp(2.*alpha.*dt*i);
%      Eyy=-kx./beta+(kx./beta).*exp(2.*beta.*dt*i);
     
%     
    
    
    Dxx=Sxx(:,i);
    Dyy=Syy(:,i);
    Dxy=Sxy(:,i);
    Dyx=Syx(:,i);
    Mdet=Dxx.*Dyy-Dxy.*Dyx;
    %if(i~=2)
    N11=Dyy./Mdet;
    N12=-Dxy./Mdet;
    N21=-Dyx./Mdet;
    N22=Dxx./Mdet;
    %else
    %   N11=N11+Dyy./Mdet; 
    %   N12=N12-Dxy./Mdet;
    %N21=N21-Dyx./Mdet;
    %N22=N22+Dxx./Mdet;
    %end
    P=N11.*XLogan.^2+N12.*XLogan.*YLogan+N21.*YLogan.*XLogan+N22.*YLogan.^2;
    
    RealP=((RealP)/(2*pi)).*Mdet.^(-1/2).*exp(-.5.*P);
%     ind=find(isnan(RealP));
%     
%     RealP(ind)=0;
%    ind=find(RealP>100000);
%    RealP(ind)=100000;
    %Ptot=RealP;
    
    
    figure(1);
    subplot(1,4,1)
    pcolor(dle); shading interp
    title('DLE')
%     subplot(1,7,2)
%     pcolor(reshape(mx,256,128)');shading interp
%     subplot(1,7,3)
%     pcolor(reshape(my,256,128)');shading interp
    subplot(1,4,2)
    pcolor(reshape(Ptot,256,256)');shading interp
    title('Sxx')
    subplot(1,4,3)
%     pcolor(reshape(sxy,256,128)');shading interp
%     subplot(1,5,4)
    pcolor(reshape(P,256,256)');shading interp
    title('P')
%     subplot(1,4,4)
%     pcolor(reshape(Exx,256,256)');shading interp
%     title('Syy')
    subplot(1,4,4)
    pcolor(reshape(RealP,256,256)');shading interp
    title('RealP')
    
    
    
    
    
    drawnow
    
    
    
    
    
    
    
    
    
    
end
    