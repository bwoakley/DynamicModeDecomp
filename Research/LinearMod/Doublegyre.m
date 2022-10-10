function [dxdt]=Doublegyre(tt,x,M,N)
global tprev;
x=reshape(x,M,N,2);
xx=x(:,:,1);
yy=x(:,:,2);
dxdt(:,:,1)=-sin(xx).*cos(yy);
dxdt(:,:,2)=cos(xx).*sin(yy);
u=dxdt(:,:,1);v=dxdt(:,:,2);dxdt=dxdt(:);
tprev=tt;
end
