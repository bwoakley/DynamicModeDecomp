close all;clear;  clc;

% figure;set(gcf,'position',[100 100 800 800])
format long; k=256;
dt=.002;Ri=1;
x=linspace(0,2*pi,k+1);x=x(1:end-1);y=x;
x1=linspace(0,2*pi,k+1);x1=x1(1:end-1);y1=x1;
[Y X]=meshgrid(y,x);
[Y1 X1]=meshgrid(y1,x1);
ss='Diff';
xx=linspace(-2*pi,4*pi,3*k+1);xx=xx(1:end-1);yy=xx;
for ri=1:Ri
xnow(:,:,ri)=X1;ynow(:,:,ri)=Y1;
end
for i=1000:1 :11000
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
    u=data(:,:,1);
    v=data(:,:,2);
    om=data(:,:,3);
    uu1=[u u u;u u u;u u u];vv1=[v v v; v v v; v v v];
    if i+1<10
    st=strcat('../Cases/',ss,'/bin0000',num2str(i+1));
    elseif i+1<100
    st=strcat('../Cases/',ss,'/bin000',num2str(i+1));
    elseif i+1<1000
    st=strcat('../Cases/',ss,'/bin00',num2str(i+1));
    elseif i+1<10000
    st=strcat('../Cases/',ss,'/bin0',num2str(i+1));
    else
    st=strcat('../Cases/',ss,'/bin',num2str(i+1));
    end
    fid=fopen(st,'rb');
    data=fread(fid,[1 1],'*float');
    data=fread(fid,[1 inf],'*double');
    fclose(fid);
    data=reshape(data,k,k,3);
    u=data(:,:,1);
    v=data(:,:,2);
    om=data(:,:,3);
    uu2=[u u u;u u u;u u u];vv2=[v v v; v v v; v v v];
    
    
% subplot(2,2,1);pcolor(x,x,u');shading interp;daspect([1 1 1]);colorbar
% subplot(2,2,2);pcolor(x,x,v');shading interp;daspect([1 1 1]);colorbar
% subplot(2,2,3);pcolor(x,x,om');shading interp;daspect([1 1 1]);colorbar
% 
% % Load saved trajectory data
%     if i<10
%     st=strcat('../Cases/',ss,'/pos0000',num2str(i));
%     elseif i<100
%     st=strcat('../Cases/',ss,'/pos000',num2str(i));
%     elseif i<1000
%     st=strcat('../Cases/',ss,'/pos00',num2str(i));
%     elseif i<10000
%     st=strcat('../Cases/',ss,'/pos0',num2str(i));
%     else
%     st=strcat('../Cases/',ss,'/pos',num2str(i));
%     end
% fid=fopen(st,'rb');
% data=fread(fid,[1 1],'*float');
% data=fread(fid,[1 inf],'*double');
% fclose(fid);
% data=reshape(data,k,k,2);
% xnow=data(:,:,1);
% ynow=data(:,:,2);
% plot(xnow,ynow,'.k');drawnow;
% Generate trajectory data

xnow=mod(xnow,2*pi);ynow=mod(ynow,2*pi);

if (mod(i,50)==0)
[uy ux]=gradient(u,y,x);
[vy vx]=gradient(v,y,x);
OM=(vx-uy);
s1=ux-vy;s2=vx+uy;
Q=(s1.^2+s2.^2-OM.^2);QQ=[Q Q Q;Q Q Q;Q Q Q];
Q0=std(Q(:));
end
for ri=1:Ri
    if (mod(i,50)==0)
Qnow(:,:,ri)=interp2(yy,xx,QQ,ynow(:,:,ri),xnow(:,:,ri));
Qa(:,:,ri,(i-1000)/50+1)=Qnow(:,:,ri)/Q0;
    end
urhs1=interp2(yy,xx,uu1,ynow(:,:,ri),xnow(:,:,ri));
vrhs1=interp2(yy,xx,vv1,ynow(:,:,ri),xnow(:,:,ri));
x1=xnow(:,:,ri)+dt*urhs1;
y1=ynow(:,:,ri)+dt*vrhs1;
xnow(:,:,ri)=x1;ynow(:,:,ri)=y1;
plot(xnow,ynow,'.k');
drawnow;
end



% xnow=x1;ynow=y1;

end
% save Qhistory.mat Qhis
% Qhist