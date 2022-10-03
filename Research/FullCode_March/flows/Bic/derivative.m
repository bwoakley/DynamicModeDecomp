function [u,v]=derivative(tt,xx,yy)
beta=.6144;
        e1=0.002;
        e2=0.3;
        Delta=sqrt(1-beta*3/2);
        c1=1/3*(1+Delta);
        c2=1/3*(1-Delta);
        k1=sqrt(6*c1);
        k2=sqrt(6*c2); 
        Om=k1*(c1-c2);
u=sech(yy).^2.*(1+2*e2*cos(k2*xx).*tanh(yy)+2*e1*cos(k1*xx-Om*tt).*tanh(yy))-c2;
v=sech(yy).^2.*(-k2*e2*sin(k2*xx)-k1*e1*sin(k1*xx-Om*tt));
end



% function [u,v]=derivative(tt,xx,yy,param)
% 
% u=sech(yy).^2.*(1+2*param.e2*cos(param.k2*xx).*tanh(yy)+2*param.e1*cos(param.k1*xx-param.Om*tt).*tanh(yy))-param.c2;
% v=sech(yy).^2.*(-param.k2*param.e2*sin(param.k2*xx)-param.k1*param.e1*sin(param.k1*xx-param.Om*tt));
% 
% end
