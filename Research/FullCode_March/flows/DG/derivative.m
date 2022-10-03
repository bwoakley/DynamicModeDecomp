function [u,v]=derivative(tt,xx,yy)

e1=3;
e2=0;%.1;
g=e1*sin(4*pi*tt)+e2*sin(2*tt);
u=-pi*sin(pi*(xx-g))*cos(pi*yy);
v=pi*cos(pi*(xx-g))*sin(pi*yy);
end



% function [u,v]=derivative(tt,xx,yy,param)
% 
% u=sech(yy).^2.*(1+2*param.e2*cos(param.k2*xx).*tanh(yy)+2*param.e1*cos(param.k1*xx-param.Om*tt).*tanh(yy))-param.c2;
% v=sech(yy).^2.*(-param.k2*param.e2*sin(param.k2*xx)-param.k1*param.e1*sin(param.k1*xx-param.Om*tt));
% 
% end
