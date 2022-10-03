% close all;clear;  clc;

figure;set(gcf,'position',[100 100 1250 550])
format long; k=256;
x=linspace(0,2,k+1);x=x(1:end-1);y=x;
[Y X]=meshgrid(y,x);
xx=[x-2*pi x x+2*pi];
%ss='FRNT/Eddy';
ss='NFlo';


t1=1001;

st=strcat('../Cases/',ss,'/bin0',num2str(t1));
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,k,k,3);
omold=data(:,:,3);
uold=data(:,:,1);
vold=data(:,:,2);

st=strcat('../Cases/',ss,'/Th0',num2str(t1));
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,k,k,1);
Cold=data(:,:,1);

dt=0.02;
for i=t1+1:1 :t1+200
    i
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
om=data(:,:,3);
u=data(:,:,1);
v=data(:,:,2);
subplot(1,2,1)
pcolor(x,x,u');shading interp;colorbar
set(gca,'fontsize',18)
daspect([1,1,1])
% caxis([.6 1.9])
title('C')
hold on;

[u(114,212) v(24,188)]


U1=[uold uold uold;
    uold uold uold;
    uold uold uold];
V1=[vold vold vold;
    vold vold vold;
    vold vold vold];

U2=[u u u;
    u u u;
    u u u];
V2=[v v v;
    v v v;
    v v v];
Uav=U1/2+U2/2;
Vav=V1/2+V2/2;

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
CC=data;
subplot(1,2,2)
pcolor(x,x,CC');shading interp;colorbar
set(gca,'fontsize',18)
daspect([1,1,1])
% caxis([.6 1.9])
title('C')
hold on;

C=contourc(x,x,CC',[.5 .5]);
ll=size(C,2);
li=0;
ld=1;
while ld<ll
    li=li+1;
L(li)=C(2,ld);
xc=C(1,ld+1:ld+L(li));
yc=C(2,ld+1:ld+L(li));
plot(xc,yc,'k')
k1x=interp2(xx,xx,U1,yc,xc);
k1y=interp2(xx,xx,V1,yc,xc);
xc1=xc+dt*k1x/2;
yc1=yc+dt*k1y/2;
k2x=interp2(xx,xx,Uav,yc1,xc1);
k2y=interp2(xx,xx,Vav,yc1,xc1);
xc2=xc+dt*k2x/2;
yc2=yc+dt*k2y/2;
k3x=interp2(xx,xx,Uav,yc2,xc2);
k3y=interp2(xx,xx,Vav,yc2,xc2);
xc3=xc+dt*k3x;
yc3=yc+dt*k3y;
k4x=interp2(xx,xx,U2,yc3,xc3);
k4y=interp2(xx,xx,V2,yc3,xc3);
xcc=xc+dt/6*(k1x+2*k2x+2*k3x+k4x);
ycc=yc+dt/6*(k1y+2*k2y+2*k3y+k4y);



plot(xcc,ycc,'r')

ld=ld+L(li)+1;
end

if i>1000
C=contourc(x,x,Cold(:,:,1)',[.5 .5]);
ll=size(C,2);
li=0;
ld=1;
while ld<ll
    li=li+1;
L(li)=C(2,ld);
xc=C(1,ld+1:ld+L(li));
yc=C(2,ld+1:ld+L(li));

plot(xc,yc,'w')

ld=ld+L(li)+1;
end
end

set(gca,'fontsize',18,'position',[.54 .05 .38 .9])
daspect([1,1,1])
hold off


subplot(1,2,1)
pcolor(x,x,om');shading interp;colorbar
daspect([1,1,1]);axis tight
set(gca,'fontsize',18,'position',[.05 .05 .38 .9])
title('\omega')
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
plot(xc,yc,'k')
k1x=interp2(xx,xx,U1,yc,xc);
k1y=interp2(xx,xx,V1,yc,xc);
xc1=xc+dt*k1x/2;
yc1=yc+dt*k1y/2;
k2x=interp2(xx,xx,Uav,yc1,xc1);
k2y=interp2(xx,xx,Vav,yc1,xc1);
xc2=xc+dt*k2x/2;
yc2=yc+dt*k2y/2;
k3x=interp2(xx,xx,Uav,yc2,xc2);
k3y=interp2(xx,xx,Vav,yc2,xc2);
xc3=xc+dt*k3x;
yc3=yc+dt*k3y;
k4x=interp2(xx,xx,U2,yc3,xc3);
k4y=interp2(xx,xx,V2,yc3,xc3);
xcc=xc+dt/6*(k1x+2*k2x+2*k3x+k4x);
ycc=yc+dt/6*(k1y+2*k2y+2*k3y+k4y);



plot(xcc,ycc,'r')

ld=ld+L(li)+1;
end

if i>1000
C=contourc(x,x,Cold(:,:,1)',[.5 .5]);
ll=size(C,2);
li=0;
ld=1;
while ld<ll
    li=li+1;
L(li)=C(2,ld);
xc=C(1,ld+1:ld+L(li));
yc=C(2,ld+1:ld+L(li));

plot(xc,yc,'w')

ld=ld+L(li)+1;
end
end
hold off

drawnow


Cold=CC;
uold=u;
vold=v;
end
