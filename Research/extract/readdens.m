close all;clear;  clc;

figure;set(gcf,'position',[100 100 800 700])
format long; k=256;
x=linspace(0,2,k+1);x=x(1:end-1);y=x;
[Y X]=meshgrid(y,x);
% for i=1000:1:1000
%     if i<10
%     st=strcat('../Cases/Diff/bin0000',num2str(i));
%     elseif i<100
%     st=strcat('../Cases/Diff/bin000',num2str(i));
%     elseif i<1000
%     st=strcat('../Cases/Diff/bin00',num2str(i));
%     else
%     st=strcat('../Cases/Diff/bin0',num2str(i));
%     end
% fid=fopen(st,'rb');
% % data=fread(fid,[1 1],'*float');
% data=fread(fid,[1 inf],'*double');
% fclose(fid);
% data=reshape(data,k,k,3);
% subplot(2,2,1)
% pcolor(x,x,data(:,:,3)');shading interp;colorbar
% daspect([1,1,1])
% set(gca,'fontsize',18)
% title('\omega')
% drawnow
% end

ss='Line';

for i=1000:1 :1010
    
    if i<10
    st=strcat('../Cases/',ss,'/bin0000',num2str(i));
    elseif i<100
    st=strcat('../Cases/',ss,'/bin000',num2str(i));
    elseif i<1000
    st=strcat('../Cases/',ss,'/bin00',num2str(i));
    elseif i<10000
    st=strcat('../Cases/',ss,'/bin0',num2str(i));
    else
    st=strcat('../Cases/',ss,'/bin',num2str(i));
    end
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);

data=reshape(data,k,k,3);

%om=data(:,:,3);
%ooo(i+1)=sum(om(:).^2);

subplot(2,2,1)
pcolor(data(:,:,3)');shading interp;
daspect([1,1,1]);axis tight
set(gca,'fontsize',18)
title('\omega')
colorbar

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
subplot(2,2,2)
pcolor(x,x,data(:,:,1)');shading interp;colorbar
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
plot(xc,yc,'k')
ld=ld+L(li)+1;

end

drawnow
hold off;

% N=data(:,:,1);B=data(:,:,2);
% BC(i-999)=sum(N(:).*B(:))/length(N(:));
% NB(i-999)=sum(N(:).*1)/length(N(:));

% N(i+1-1000)=data(316,198,1);
% P(i+1-1000)=data(316,198,2);
% Z(i+1-1000)=data(316,198,3);
end
% t=0:.2:200*.2;
% % subplot(2,2,1)
% % plot(t,BC,t,NB)
% Adv=sum(BC-NB)*.2
