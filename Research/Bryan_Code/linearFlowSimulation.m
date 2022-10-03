close all;clear;  clc;

figure;set(gcf,'position',[100 100 800 700])
format long; k=256;
x=linspace(0,2,k+1);x=x(1:end-1);y=x;
[Y X]=meshgrid(y,x);


ss='Line';

v = VideoWriter('newfile.avi');
open(v)

for i=1000:1 :1050
     
    if i<10
    st=strcat('../Cases/',ss,'/Th0000',num2str(i));
    elseif i<100
    st=strcat('../Cases/',ss,'/Th000',num2str(i));
    elseif i<1000
    st=strcat('../Cases/',ss,'/Th00',num2str(i));
    elseif i<10000
    st=strcat('../Cases/',ss,'/Th0',num2str(i));
    else
    st=strcat('../Cases/',ss,'/Th',num2str(i));
    end


fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,k,k,1);
% subplot(2,2,2)
%pcolor(x,x,data(:,:,1)');shading interp;colorbar
set(gca,'fontsize',18)
daspect([1,1,1])
% caxis([.6 1.9])
title('C')
hold on;
C=contourc(x,x,data(:,:,1)',[.5 .5]);

ll=size(C,2);
li=0;
ld=1;
while ld<ll
    li=li+1;
    L(li)=C(2,ld);
    xc=C(1,ld+1:ld+L(li));
    yc=C(2,ld+1:ld+L(li));
    for ki=1:length(xc)-1
        xkc(ki*2-1)=xc(ki);xkc(ki*2)=(xc(ki)+xc(ki+1))/2;
        ykc(ki*2-1)=yc(ki);ykc(ki*2)=(yc(ki)+yc(ki+1))/2;
    end
    %plot(xc,yc,'k');
    ld=ld+L(li)+1;
end

%drawnow
%Instead of drawnow, save the movie in F
%F(i) = getframe;
%movie2avi(F,'myavifile.avi','Compression','Cinepak')

writeVideo(v,data(:,:,1)')


hold off;


end

close(v)
