% close all;clear;  clc;

figure;set(gcf,'position',[100 100 800 350])
format long; k=64;
x=linspace(0,2*pi,k+1);x=x(1:end-1);y=x;
[Y X]=meshgrid(y,x);

st='../Cases/BSTA/Da4Th';
fid=fopen(st,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);
data=reshape(data,k,k,251);

for i=1:251
pcolor(data(:,:,i)');shading interp;colorbar

dd=data(:,:,i);
daspect([1,1,1])
drawnow
%gap(i)=max(dd(:))-min(dd(:));
end
%plot(gap)
