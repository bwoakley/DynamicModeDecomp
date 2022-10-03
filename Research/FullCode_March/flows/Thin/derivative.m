function [u,v]=derivative(tt,xx,yy)

u=.0001.*(yy-2);
v=0.*xx;
end



% function [u,v]=derivative(tt,xx,yy,param)
% 
% u=sech(yy).^2.*(1+2*param.e2*cos(param.k2*xx).*tanh(yy)+2*param.e1*cos(param.k1*xx-param.Om*tt).*tanh(yy))-param.c2;
% v=sech(yy).^2.*(-param.k2*param.e2*sin(param.k2*xx)-param.k1*param.e1*sin(param.k1*xx-param.Om*tt));
% 
% end
