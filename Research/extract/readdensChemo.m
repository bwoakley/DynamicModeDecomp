close all;clear;  clc;

figure;set(gcf,'position',[100 100 800 350])
format long; k=256;
x=linspace(0,2*pi,k+1);x=x(1:end-1);y=x;
[Y X]=meshgrid(y,x);
% for i=1000:1:1010
%     if i<10
%     st=strcat('../Cases/diff/bin000',num2str(i));
%     elseif i<100
%     st=strcat('../Cases/diff/bin00',num2str(i));
%     elseif i<1000
%     st=strcat('../Cases/diff/bin0',num2str(i));
%     else
%     st=strcat('../Cases/diff/bin',num2str(i));
%     end
% fid=fopen(st,'rb');
% data=fread(fid,[1 1],'*float');
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

ss='BSTA/Da8/Hyp';
ss1='BSTA/Da8/Eddy';

for i=1200:1 :1200
    i
%     if i<10
%     st=strcat('../Cases/',ss1,'/bin0000',num2str(i));
%     elseif i<100
%     st=strcat('../Cases/',ss1,'/bin000',num2str(i));
%     elseif i<1000
%     st=strcat('../Cases/',ss1,'/bin00',num2str(i));
%     elseif i<10000
%     st=strcat('../Cases/',ss1,'/bin0',num2str(i));
%     else
%     st=strcat('../Cases/',ss1,'/bin',num2str(i));
%     end
% fid=fopen(st,'rb');
% data=fread(fid,[1 1],'*float');
% data=fread(fid,[1 inf],'*double');
% fclose(fid);
% data=reshape(data,k,k,3);
% om=data(:,:,3);
% ooo(i+1)=sum(om(:).^2);
% 
% subplot(2,2,1)
% pcolor(data(:,:,3)');shading interp;
% daspect([1,1,1]);axis tight
% set(gca,'fontsize',18)
% title('\omega')
% colorbar
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
data=reshape(data,k,k,3);
subplot(1,2,1)
pcolor(x,x,data(:,:,1)');shading interp;colorbar
set(gca,'fontsize',18)
set(gca,'xtick',[0 2 4 6]);
set(gca,'ytick',[0 2 4 6]);
daspect([1,1,1])
% caxis([.6 1.9])
title('C')
xlabel('X')
ylabel('Y')
    if i<10
    st=strcat('../Cases/',ss1,'/Th0000',num2str(i));
    elseif i<100
    st=strcat('../Cases/',ss1,'/Th000',num2str(i));
    elseif i<1000
    st=strcat('../Cases/',ss1,'/Th00',num2str(i));
    elseif i<10000
    st=strcat('../Cases/',ss1,'/Th0',num2str(i));
    else
    st=strcat('../Cases/',ss1,'/Th',num2str(i));
    end
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,k,k,3);subplot(1,2,2)
pcolor(x,x,data(:,:,1)');shading interp;colorbar
daspect([1,1,1])
set(gca,'fontsize',18)
set(gca,'xtick',[0 2 4 6]);
set(gca,'ytick',[0 2 4 6]);
% caxis([.3 2.1])
title('C')
xlabel('X')
ylabel('Y')
 % % subplot(2,2,4)
% % pcolor(x,x,data(:,:,3)');shading interp;colorbar
% % daspect([1,1,1])
% % set(gca,'fontsize',18)
% % caxis([.0 .5])
% % title('Z')
drawnow

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