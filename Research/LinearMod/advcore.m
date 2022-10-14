clear all;close all;
clc;
format long;

global tprev lambda1 mu1 LAM MU;
global thm thp tau_e phitot maxphi phi_t phim_t phip_t T Tflag Period;
global tau_h sigma x_t y_t phi0o dir start_time maxtau_e tau_eseg;

%Parameters for Bickley jet
%Full period
Period=1

%Case 
CaseSTART=0; %Start at t=0, find LCS at different base times
CaseEND=0;   %End at t=1 period.
dir=1; %%Direction of integration, +1 F, -1 B

BASESTEP=0.2*Period; %basetime interval
INTTIME=.1*Period; %integration time

%Discretization
x=linspace(0,2*pi,256+1);
y=linspace(0,1*pi,128+1);

% x=x(150);y=y(80);
% th=(0:2:359)/360*2*pi;
% r=1e-3:1e-3:1e-1;
% [TH R]=meshgrid(th,r);
% xr=x+R.*cos(TH);
% yr=y+R.*sin(TH);



% x=-1:.1:1;y=0:.1:2;
%x=2.62;y=-.085;
M=length(x);N=length(y);
[YST XST]=meshgrid(y,x); 
%XST=xr;
%YST=yr;
Diff=CaseEND-CaseSTART;
N_Basetime=Diff/(BASESTEP/Period)+1 %%# of base times for animation
%%Set integration time
start_time=CaseSTART;
xt=XST;
yt=YST;
for ti=1:N_Basetime
    Current_Base=(ti-1)*BASESTEP/Period
    
    xnow = XST;%-1e-6;
    ynow = YST;%+1e-6;
    tj=0;
             
    Tflag=0;

    
%Solve xdot=u with ode45 to set accuracy
%Have to do this iteratively for smaller intervals for space
    for timeunit=1:INTTIME/Period*20
        current_unit=timeunit/20
        Tspan=([0 Period/20]+Tflag*Period/20)*dir+start_time
        IC(:,:,1)=xnow;
        IC(:,:,2)=ynow;
        options = odeset('RelTol',1e-14,'AbsTol',1e-14,'Maxstep',0.01);
        [TT, xpos]=ode45(@(t,xpos) Doublegyre(t,xpos,M,N),...
            Tspan,IC);%,options);
        xfinal=reshape(xpos(end,:),M,N,2);
        xnow=xfinal(:,:,1);
        ynow=xfinal(:,:,2);
        Tflag=Tflag+1;
        xt(:,:,timeunit+1)=xnow;
        yt(:,:,timeunit+1)=ynow;
        
%         if mod(timeunit,INTTIME/Period*20)==0
        %xa=xnow;ya=ynow;
        xa = x; ya = y;
        dle_2d;
%         end
    end
%Renew base time    
    start_time=start_time+BASESTEP;
end
