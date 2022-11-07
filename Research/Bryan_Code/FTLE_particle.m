clear all;
close all;
clc;
format long;

kk = 256;

%Length of domain
L = 3;
%Full period or time unit
Period=1;
%Find a point of high, med, low stretching. Here are the coords:
rowIndexHigh = 201; %y-coord index of High
colIndexHigh = 73; %x-coord index of High
rowIndexMed = 120;
colIndexMed = 143;
rowIndexLow = 65;
colIndexLow = 137;

%Case 
CaseSTART=200-.2; %Start at t=0, find LCS at different base times
CaseEND=200-.2;   %End at t=1 period.
dir=1; %%Direction of integration, +1 F, -1 B

timestepsize=.02; 
deltat = dir*timestepsize; %%integration step size with direction

BASESTEP=0.2*Period; %basetime interval
%INTTIME=2*Period; %integration time. Right now, this goes to time t = 2. Let's reduce it:
INTTIME=49*abs(deltat);
%INTTIME=.99*Period;


%Discretization
x=linspace(0,L,256*3+1); 
y=linspace(0,L,256*3+1);
xx=x(1:end-1);
%Above is a large grid, but my data is not defined on a 768 grid... 
%Instead we will define it on a 256 grid that still spans the whole length L
xcoarse = linspace(0,L,256+1);
ycoarse = linspace(0,L,256+1);
xxcoarse = xcoarse(1:end-1);
yycoarse = ycoarse(1:end-1);
[XSTcoarse, YSTcoarse]=meshgrid(xxcoarse,yycoarse); 

inter='*cubic';
%inter='*linear';
%inter='linear';

xa=x(257:512);ya=y(257:512);
[XST, YST]=meshgrid(xa,ya); 

Diff=CaseEND-CaseSTART;
N_Basetime=Diff/(BASESTEP/Period)+1; %%# of base times for animation
%%Set integration time
start_time=CaseSTART;

no_of_steps = INTTIME/abs(deltat);
ft=.2;
for ti=1:N_Basetime
    Current_Base=(ti-1)*BASESTEP/Period;
    
    xnow = XST;%-1e-6;
    ynow = YST;%+1e-6;

    %Keep track of the position of the high, med, low particles.
    xCoordHistoryHigh = zeros(1,no_of_steps);
    yCoordHistoryHigh = zeros(1,no_of_steps);
    xCoordHistoryMed = zeros(1,no_of_steps);
    yCoordHistoryMed = zeros(1,no_of_steps);
    xCoordHistoryLow = zeros(1,no_of_steps);
    yCoordHistoryLow = zeros(1,no_of_steps);

    %Keep track of the FTLE for high, med, low
    FTLEHistoryHigh = zeros(1,no_of_steps);
    FTLEHistoryMed = zeros(1,no_of_steps);
    FTLEHistoryLow = zeros(1,no_of_steps);


    for i = 1:no_of_steps
        i;
        current_time = start_time + (i-1)*deltat;   
        x0 = xnow;
        y0 = ynow;
        t0 = current_time;


        %Update the current particle location
        xCoordHistoryHigh(i) = xnow(rowIndexHigh,colIndexHigh);
        yCoordHistoryHigh(i) = ynow(rowIndexHigh,colIndexHigh);
        xCoordHistoryMed(i) = xnow(rowIndexMed,colIndexMed);
        yCoordHistoryMed(i) = ynow(rowIndexMed,colIndexMed);
        xCoordHistoryLow(i) = xnow(rowIndexLow,colIndexLow);
        yCoordHistoryLow(i) = ynow(rowIndexLow,colIndexLow);


        %data file name
        ss='Turb';

        ii = i + 1000;
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

        
       

        xm=mod(x0,L);ym=mod(y0,L);
        %urhs1 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        %vrhs1 = interp2(xx,xx,vinterp,ym,xm,inter,0);          
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
        urhs1 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
        vrhs1 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);


%*******************************************************************
        x1 = x0 + urhs1*deltat/2;
        y1 = y0 + vrhs1*deltat/2;
        %min(min(x1)) %Yes, x1 does leave the set [1,2]^2 a bit. So, we should mod it back.
        %max(max(x1))
        t1 = t0 + deltat/2;
        index1n = floor(t1/ft) + 1;
        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %checkload;
        end 
        uinterp = u1 + (t1 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t1 - (index1-1)*ft)*(v2-v1)/ft;

        %xm=mod(x0,L);ym=mod(y0,L);
        xm=mod(x1,L);ym=mod(y1,L);   %***********xm should update to mod(x1,L)   
        
        %urhs2 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        %vrhs2 = interp2(xx,xx,vinterp,ym,xm,inter,0);
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
        urhs2 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
        vrhs2 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %*************************************************************
        x2 = x0 + urhs2*deltat/2;
        y2 = y0 + vrhs2*deltat/2;
        t2 = t0 + deltat/2;   

        %xm=mod(x0,L);ym=mod(y0,L);
        xm=mod(x2,L);ym=mod(y2,L); 

        %urhs3 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        %vrhs3 = interp2(xx,xx,vinterp,ym,xm,inter,0);   
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
        urhs3 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
        vrhs3 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
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
        
        %xm=mod(x0,L);ym=mod(y0,L);
        xm=mod(x3,L);ym=mod(y3,L);  

        %urhs4 = interp2(xx,xx,uinterp,ym,xm,inter,0);
        %vrhs4 = interp2(xx,xx,vinterp,ym,xm,inter,0); 
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
        urhs4 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
        vrhs4 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %********************************************************************   
    
        dx = (deltat/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
        dy = (deltat/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
  
       %First call FTLE code to use xbefore = xnow and xafter
       xbefore = xnow;
       xafter = xnow + dx;
       ybefore = ynow;
       yafter = ynow + dy;
       
       dle_2d;

       %Then update xnow:
        xnow = xnow + dx;
        ynow = ynow + dy;
        
       
       
%         figure;
%         plot(xnow,ynow,'.k'); daspect([1 1 1]); 
%         drawnow;
%         pause(.1)
        
        %Now compute the FTLE at every time step

    end

    figure;
    c = linspace(1,10,length(xCoordHistoryHigh));
    scatter(xCoordHistoryHigh,yCoordHistoryHigh,[],c,'filled')
    hold on;
    scatter(xCoordHistoryMed,yCoordHistoryMed,[],c,'d')
    hold on;
    scatter(xCoordHistoryLow,yCoordHistoryLow,[],c)
    hold off;
    xlim([1,2])
    ylim([1,2])
    title('particleHistory. solid = high, diamond = med, open circle = low')

    figure;
    subplot(3,1,1)
    plot(FTLEHistoryHigh)
    title('FTLE High')
    subplot(3,1,2)
    plot(FTLEHistoryMed)
    title('FTLE Med')
    subplot(3,1,3)
    plot(FTLEHistoryLow)
    title('FTLE Low')

    
%Renew base time    
    start_time=start_time+BASESTEP;
end

%ss=strcat('save FTLEBTurb.mat FTLE');
%eval(ss);