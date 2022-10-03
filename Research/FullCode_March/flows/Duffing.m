function [dxdt]=Duffing(t,x,M,N)

e1=.1;

x=reshape(x,M,N,2);
xx=x(:,:,1);
yy=x(:,:,2);
dxdt(:,:,1)=yy;
dxdt(:,:,2)=xx-xx.^3+.1*sin(t);
dxdt=dxdt(:);
end
