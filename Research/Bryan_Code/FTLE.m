clear all;close all;
clc;
format long;

global tprev lambda1 mu1 LAM MU;
global thm thp tau_e phitot maxphi phi_t phim_t phip_t T Tflag Period;
global tau_h sigma x_t y_t phi0o dir start_time maxtau_e tau_eseg;

kk = 256;

%Parameters for Bickley jet
%Full period or time unit
Period=1

%Case 
CaseSTART=200-.2; %Start at t=0, find LCS at different base times
CaseEND=200-.2;   %End at t=1 period.
dir=1; %%Direction of integration, +1 F, -1 B

BASESTEP=0.2*Period; %basetime interval
INTTIME=2*Period; %integration time
timestepsize=.02; 
deltat = dir*timestepsize; %%integration step size with direction

%Discretization
x=linspace(0,6*pi,256*3+1);
y=linspace(0,6*pi,256*3+1);
xx=x(1:end-1);
% x=x(150);y=y(80);
% th=(0:2:359)/360*2*pi;
% r=1e-3:1e-3:1e-1;
% [TH R]=meshgrid(th,r);
% xr=x+R.*cos(TH);
% yr=y+R.*sin(TH);

inter='*cubic';

% x=-1:.1:1;y=0:.1:2;
%x=2.62;y=-.085;
M=length(x);N=length(y);
xa=x(257:512);ya=y(257:512);
[YST XST]=meshgrid(y(257:512),x(257:512)); 
%XST=xr;
%YST=yr;
Diff=CaseEND-CaseSTART;
N_Basetime=Diff/(BASESTEP/Period)+1 %%# of base times for animation
%%Set integration time
start_time=CaseSTART;
xt=XST;
yt=YST;
no_of_steps = INTTIME/abs(deltat);
ft=.2;
for ti=1:N_Basetime
    Current_Base=(ti-1)*BASESTEP/Period
    
    xnow = XST;%-1e-6;
    ynow = YST;%+1e-6;

    for i = 1:no_of_steps
        i
        current_time = start_time + (i-1)*deltat;   
        x0 = xnow;
        y0 = ynow;
        t0 = current_time;

        %data file name
        ss='Line';
        ii = i + 999;
        st=strcat('../Cases/',ss,'/bin0',num2str(ii));
        fid=fopen(st,'rb');
        data=fread(fid,[1 1],'*float');
        data=fread(fid,[1 inf],'*double');
        fclose(fid);
        data=reshape(data,kk,kk,3);
        u1 = data(:,:,1)' ;
        v1 = data(:,:,2)' ;

        ii2 = ii + 1;
        st=strcat('../Cases/',ss,'/bin0',num2str(ii2));
        fid=fopen(st,'rb');
        data=fread(fid,[1 1],'*float');
        data=fread(fid,[1 inf],'*double');
        fclose(fid);
        data=reshape(data,kk,kk,3);
        u2 = data(:,:,1)' ;
        v2 = data(:,:,2)' ;

        loaded = [-1000 -1000];
        index1=-10000;
        index1n = floor(t0/ft) + 1;

        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %run this check to save loading time
            %checkload;
        end     
        uinterp = u1 + (t0 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t0 - (index1-1)*ft)*(v2-v1)/ft;

        xm=mod(x0,6*pi);ym=mod(y0,6*pi);
        urhs1 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        vrhs1 = interp2(xx,xx,vinterp,ym,xm,inter,0);
%*******************************************************************
        x1 = x0 + urhs1*deltat/2;
        y1 = y0 + vrhs1*deltat/2;
        t1 = t0 + deltat/2;
        index1n = floor(t1/ft) + 1;
        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %checkload;
        end 
        uinterp = u1 + (t1 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t1 - (index1-1)*ft)*(v2-v1)/ft;
        xm=mod(x0,6*pi);ym=mod(y0,6*pi);        
        urhs2 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        vrhs2 = interp2(xx,xx,vinterp,ym,xm,inter,0);
 %*************************************************************
        x2 = x0 + urhs2*deltat/2;
        y2 = y0 + vrhs2*deltat/2;
        t2 = t0 + deltat/2;   
        xm=mod(x0,6*pi);ym=mod(y0,6*pi);        
        urhs3 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        vrhs3 = interp2(xx,xx,vinterp,ym,xm,inter,0);   
 %****************************************************************
        x3 = x0 + urhs3*deltat;
        y3 = y0 + vrhs3*deltat;
        t3 = t0 + deltat;
        index1n = floor(t3/ft) + 1;
        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %checkload;
        end 
        
        uinterp = u1 + (t3 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t3 - (index1-1)*ft)*(v2-v1)/ft;
        xm=mod(x0,6*pi);ym=mod(y0,6*pi);        
        urhs4 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        vrhs4 = interp2(xx,xx,vinterp,ym,xm,inter,0);  
 %********************************************************************   
    
        dx = (deltat/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
        dy = (deltat/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
        xnow = xnow + dx;
        ynow = ynow + dy;

    end
    dle_2d;
    

    
%Renew base time    
    start_time=start_time+BASESTEP;
end
ss=strcat('save FTLEBTurb.mat FTLE');
eval(ss);