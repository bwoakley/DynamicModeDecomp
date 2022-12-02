clear all;
close all;
clc;
format long;


kk = 256;

%Length of domain
L = 3;
%Full period or time unit
Period=1;
%Find a point of high, med, low stretching. Here are the coords, chosen
%after 99 time steps:
rowIndexHigh = 151; %y-coord index of High
colIndexHigh = 95; %x-coord index of High
rowIndexMed = 127;
colIndexMed = 94;
rowIndexLow = 144;
colIndexLow = 64;

%Case 
CaseSTART=200-.2; %Start at t=0, find LCS at different base times
CaseEND=200-.2;   %End at t=1 period.
dir=1; %%Direction of integration, +1 F, -1 B

timestepsize=.02; 
deltat = dir*timestepsize; %%integration step size with direction

BASESTEP=0.2*Period; %basetime interval
%INTTIME=2*Period; %integration time. Right now, this goes to time t = 2. Let's reduce it:
%INTTIME=99*abs(deltat);
INTTIME=3*abs(deltat);
%INTTIME=.99*Period;


%Discretization
x=linspace(0,L,256*3+1); 
y=linspace(0,L,256*3+1);
xx=x(1:end-1);
yy=y(1:end-1);

%Above is a large grid, but my data is not defined on a 768 grid... 

    %Instead we will define it on a 256 grid that still spans the whole length L
    xcoarse = linspace(0,L,256+1);
    ycoarse = linspace(0,L,256+1);
    xxcoarse = xcoarse(1:end-1);
    yycoarse = ycoarse(1:end-1);
    [XSTcoarse, YSTcoarse]=meshgrid(xxcoarse,yycoarse); 

    %Instead of coarse --> interp to xa, lets just take 3 copies of  the
    %flow...


inter='*cubic';
%inter='*linear';
%inter='linear';

xa=x(257:512);ya=y(257:512);
[XST, YST]=meshgrid(xa,ya); 

Diff=CaseEND-CaseSTART;
N_Basetime=Diff/(BASESTEP/Period)+1; %%# of base times for animation
%%Set integration time
start_time=CaseSTART;




no_of_steps = round(INTTIME/abs(deltat));
ft=.2;
for ti=1:N_Basetime
    Current_Base=(ti-1)*BASESTEP/Period;
    
    xnow = XST;%-1e-6;
    ynow = YST;%+1e-6;

    %Keep track of the position of the high, med, low particles.
    xCoordHistoryHigh = zeros(1,no_of_steps+1);
    yCoordHistoryHigh = zeros(1,no_of_steps+1);
    xCoordHistoryMed = zeros(1,no_of_steps+1);
    yCoordHistoryMed = zeros(1,no_of_steps+1);
    xCoordHistoryLow = zeros(1,no_of_steps+1);
    yCoordHistoryLow = zeros(1,no_of_steps+1);

    %Keep track of the FTLE for high, med, low
    FTLEHistoryHigh = zeros(1,no_of_steps);
    FTLEHistoryMed = zeros(1,no_of_steps);
    FTLEHistoryLow = zeros(1,no_of_steps);

    %Keep track of the FTLE for high, med, low. Using the 4 points
    %algorithm
    FTLEHistoryHighPoints = zeros(1,no_of_steps);
    FTLEHistoryMedPoints = zeros(1,no_of_steps);
    FTLEHistoryLowPoints = zeros(1,no_of_steps);

    %Keep track of the Strain for high, med, low
    StrainHistoryHigh = zeros(1,no_of_steps);
    StrainHistoryMed = zeros(1,no_of_steps);
    StrainHistoryLow = zeros(1,no_of_steps);

    %Update the current particle location
    xCoordHistoryHigh(1) = XST(rowIndexHigh,colIndexHigh);
    yCoordHistoryHigh(1) = YST(rowIndexHigh,colIndexHigh);
    xCoordHistoryMed(1) = XST(rowIndexMed,colIndexMed);
    yCoordHistoryMed(1) = YST(rowIndexMed,colIndexMed);
    xCoordHistoryLow(1) = XST(rowIndexLow,colIndexLow);
    yCoordHistoryLow(1) = YST(rowIndexLow,colIndexLow);

    %Keep track of the FTLE for high, med, low
    FTLENetHistoryHigh = zeros(1,no_of_steps);
    FTLENetHistoryMed = zeros(1,no_of_steps);
    FTLENetHistoryLow = zeros(1,no_of_steps);

    %Keep track of a circle of points near a chosen trajectory.
    caseSelectCircle = 3;   %1 = high, 2 = med, 3 = low
    N = 1000;               %Number of angles
    hh = 1/(256*10);        %Radius of initial circle
    vxHistory = zeros(N,no_of_steps+1);
    vyHistory = zeros(N,no_of_steps+1);
    vxHistoryPlot = zeros(N,no_of_steps+1);
    vyHistoryPlot = zeros(N,no_of_steps+1);
    thetaRange = linspace(0, 2*pi, N);
    for theta = 1:N
        vxHistoryPlot(theta,1) = cos(thetaRange(theta));
        vyHistoryPlot(theta,1) = sin(thetaRange(theta));
        if caseSelectCircle == 1
            vxHistory(theta,1) = hh*vxHistoryPlot(theta,1)+xCoordHistoryHigh(1);
            vyHistory(theta,1) = hh*vyHistoryPlot(theta,1)+yCoordHistoryHigh(1);
        elseif caseSelectCircle == 2
            vxHistory(theta,1) = hh*vxHistoryPlot(theta,1)+xCoordHistoryMed(1);
            vyHistory(theta,1) = hh*vyHistoryPlot(theta,1)+yCoordHistoryMed(1);
        else
            vxHistory(theta,1) = hh*vxHistoryPlot(theta,1)+xCoordHistoryLow(1);
            vyHistory(theta,1) = hh*vyHistoryPlot(theta,1)+yCoordHistoryLow(1);
        end
    end
    v = VideoWriter('LowStretch.avi');
    open(v)

    for i = 1:no_of_steps
        i;
        current_time = start_time + (i-1)*deltat;   
        x0 = xnow;
        y0 = ynow;
        t0 = current_time;


        


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

%         if i == 1
%             u1temp = u1;
%             v1temp = v1;
%             
%         else
%             u1 = u1temp;
%             v1 = v1temp;
%         end

%         u1 = -(XST-1.5);
%         v1 = YST-1.5;

%         figure;
%         quiver(u1,v1)

        u1Data=u1;
        v1Data=v1;

        u1Extend = [u1, u1, u1; u1, u1, u1; u1, u1, u1];
        v1Extend = [v1, v1, v1; v1, v1, v1; v1, v1, v1];


        u1 = u1Extend;
        v1 = v1Extend;




        ii2 = ii + 1;
        st=strcat('../Cases/',ss,'/bin0',num2str(ii2));
        fid=fopen(st,'rb');
        data=fread(fid,[1 1],'*float');
        data=fread(fid,[1 inf],'*double');
        fclose(fid);
        data=reshape(data,kk,kk,3);
        u2 = data(:,:,1)' ;
        v2 = data(:,:,2)' ;

%         u2 = -(XST-1.5);
%         v2 = YST-1.5;
        
%         u2 = u1temp;
%         v2 = v1temp;
%         errorFlow =  max( max ( abs(u1temp - u2) + abs(v1temp - v2)));


        u2Extend = [u2, u2, u2; u2, u2, u2; u2, u2, u2];
        v2Extend = [v2, v2, v2; v2, v2, v2; v2, v2, v2];
        u2 = u2Extend;
        v2 = v2Extend;

        loaded = [-1000 -1000];
        index1=-10000;
        index1n = floor(t0/ft) + 1;

        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %run this check to save loading time
            %checkload;
        end     
        %t0 - (index1-1)*ft
        uinterp = u1 + (t0 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t0 - (index1-1)*ft)*(v2-v1)/ft;

%         figure;
%         quiver(uinterp,vinterp)
%         title('flow')
       
        

        xm=mod(x0,L);ym=mod(y0,L);
        urhs1 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs1 = interp2(xx,yy,vinterp,xm,ym,inter,0);          
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs1 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs1 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
       

%         figure;
%         quiver(urhs1,vrhs1)
%         title('flow')
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
        
        urhs2 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs2 = interp2(xx,yy,vinterp,xm,ym,inter,0);
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs2 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs2 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %*************************************************************
        x2 = x0 + urhs2*deltat/2;
        y2 = y0 + vrhs2*deltat/2;
        t2 = t0 + deltat/2;   

        %xm=mod(x0,L);ym=mod(y0,L);
        xm=mod(x2,L);ym=mod(y2,L); 

        urhs3 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs3 = interp2(xx,yy,vinterp,xm,ym,inter,0);   
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs3 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs3 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
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

        urhs4 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs4 = interp2(xx,yy,vinterp,xm,ym,inter,0); 
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs4 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs4 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %********************************************************************   
    
        dx = (deltat/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
        dy = (deltat/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
  
        
%         if i == 20
%             figure;
%             quiver(u1,u2)
%             title('flow')
% 
%             figure;
%             quiver(dx,dy)
%             title('displacement')

       %Then update xnow:
        xnow = xnow + dx;
        ynow = ynow + dy;
        


        %Update the current particle location
        xCoordHistoryHigh(i+1) = xnow(rowIndexHigh,colIndexHigh);
        yCoordHistoryHigh(i+1) = ynow(rowIndexHigh,colIndexHigh);
        xCoordHistoryMed(i+1) = xnow(rowIndexMed,colIndexMed);
        yCoordHistoryMed(i+1) = ynow(rowIndexMed,colIndexMed);
        xCoordHistoryLow(i+1) = xnow(rowIndexLow,colIndexLow);
        yCoordHistoryLow(i+1) = ynow(rowIndexLow,colIndexLow);


        dle_2d;
    
        rescaleDle = INTTIME*dle/(deltat*i);

        FTLENetHistoryHigh(i) = rescaleDle(rowIndexHigh,colIndexHigh);
        FTLENetHistoryMed(i) = rescaleDle(rowIndexMed,colIndexMed);
        FTLENetHistoryLow(i) = rescaleDle(rowIndexLow,colIndexLow);

%*******************************************************************
%*******************************************************************


        %Now compute the FTLE at every time step
        %This requires us to take a unif grid at this time step, and
        %compute dx:


        uinterp = u1 + (t0 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t0 - (index1-1)*ft)*(v2-v1)/ft;

        
       

        xm=mod(XST,L);ym=mod(YST,L);
        urhs1 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs1 = interp2(xx,yy,vinterp,xm,ym,inter,0);          
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs1 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs1 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);


%*******************************************************************
        x1 = XST + urhs1*deltat/2;
        y1 = YST + vrhs1*deltat/2;
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

        %xm=mod(XST,L);ym=mod(YST,L);
        xm=mod(x1,L);ym=mod(y1,L);   %***********xm should update to mod(x1,L)   
        
        urhs2 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs2 = interp2(xx,yy,vinterp,xm,ym,inter,0);
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs2 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs2 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %*************************************************************
        x2 = XST + urhs2*deltat/2;
        y2 = YST + vrhs2*deltat/2;
        t2 = t0 + deltat/2;   

        %xm=mod(XST,L);ym=mod(YST,L);
        xm=mod(x2,L);ym=mod(y2,L); 

        urhs3 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs3 = interp2(xx,yy,vinterp,xm,ym,inter,0);   
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs3 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs3 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %****************************************************************
        x3 = XST + urhs3*deltat;
        y3 = YST + vrhs3*deltat;
        t3 = t0 + deltat;
        index1n = floor(t3/ft) + 1;
        if index1n~=index1
            index1 = index1n;
            index2 = index1 + 1;
            %checkload;
        end 
        
        uinterp = u1 + (t3 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t3 - (index1-1)*ft)*(v2-v1)/ft;
        
        %xm=mod(XST,L);ym=mod(YST,L);
        xm=mod(x3,L);ym=mod(y3,L);  

        urhs4 = interp2(xx,yy,uinterp,xm,ym,inter,0);
        vrhs4 = interp2(xx,yy,vinterp,xm,ym,inter,0); 
            %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
            %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                %urhs4 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                %vrhs4 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
 %********************************************************************   
    
        dx = (deltat/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
        dy = (deltat/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
       

       %First call FTLE code to use xbefore = xnow and xafter
       xbefore = XST;
       xafter = XST + dx;
       ybefore = YST;
       yafter = YST + dy;
       
       dle_2d_particle;

       FTLEHistoryHigh(i) = interp2(xa,ya,dle,xCoordHistoryHigh(i),yCoordHistoryHigh(i),inter,0);
       FTLEHistoryMed(i) = interp2(xa,ya,dle,xCoordHistoryMed(i),yCoordHistoryMed(i),inter,0);
       FTLEHistoryLow(i) =  interp2(xa,ya,dle,xCoordHistoryLow(i),yCoordHistoryLow(i),inter,0);

       
        
        %*******************************************************************
        %*******************************************************************


        %This time, compute the FTLE at every time step a different way.
        %Consider 4 points close to my trajectory.
        h = 1/(256*30);
        %Select High, Med, Low
        for caseSelect = 1:3
            if caseSelect == 1
                pNow = [xCoordHistoryHigh(i), yCoordHistoryHigh(i)];
            elseif caseSelect == 2
                pNow = [xCoordHistoryMed(i), yCoordHistoryMed(i)];
            else
                pNow = [xCoordHistoryLow(i), yCoordHistoryLow(i)];
            end
    
            pLeft = pNow+[-h,0];
            pRight = pNow+[h,0];
            pUp = pNow+[0,h];
            pDown = pNow+[0,-h];
       
            pointsX = [pLeft(1), pNow(1), pRight(1)];
            pointsY = [pDown(2), pNow(2), pUp(2)];
    
            [pointsXX, pointsYY]=meshgrid(pointsX,pointsY); 
            
    
            %Now evolve the 3 by 3 grid of points 1 time step:
    
    
            uinterp = u1 + (t0 - (index1-1)*ft)*(u2-u1)/ft;
            vinterp = v1 + (t0 - (index1-1)*ft)*(v2-v1)/ft;
    
            
           
    
            xm=mod(pointsXX,L);ym=mod(pointsYY,L);
            urhs1 = interp2(xx,yy,uinterp,xm,ym,inter,0);
            vrhs1 = interp2(xx,yy,vinterp,xm,ym,inter,0);          
                %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
                %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                    %urhs1 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                    %vrhs1 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
    
    
    %*******************************************************************
            x1 = pointsXX + urhs1*deltat/2;
            y1 = pointsYY + vrhs1*deltat/2;
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
    
            %xm=mod(pointsXX,L);ym=mod(pointsYY,L);
            xm=mod(x1,L);ym=mod(y1,L);   %***********xm should update to mod(x1,L)   
            
            urhs2 = interp2(xx,yy,uinterp,xm,ym,inter,0);
            vrhs2 = interp2(xx,yy,vinterp,xm,ym,inter,0);
                %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
                %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                    %urhs2 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                    %vrhs2 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
     %*************************************************************
            x2 = pointsXX + urhs2*deltat/2;
            y2 = pointsYY + vrhs2*deltat/2;
            t2 = t0 + deltat/2;   
    
            %xm=mod(pointsXX,L);ym=mod(pointsYY,L);
            xm=mod(x2,L);ym=mod(y2,L); 
    
            urhs3 = interp2(xx,yy,uinterp,xm,ym,inter,0);
            vrhs3 = interp2(xx,yy,vinterp,xm,ym,inter,0);   
                %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
                %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                    %urhs3 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                    %vrhs3 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
     %****************************************************************
            x3 = pointsXX + urhs3*deltat;
            y3 = pointsYY + vrhs3*deltat;
            t3 = t0 + deltat;
            index1n = floor(t3/ft) + 1;
            if index1n~=index1
                index1 = index1n;
                index2 = index1 + 1;
                %checkload;
            end 
            
            uinterp = u1 + (t3 - (index1-1)*ft)*(u2-u1)/ft;
            vinterp = v1 + (t3 - (index1-1)*ft)*(v2-v1)/ft;
            
            %xm=mod(pointsXX,L);ym=mod(pointsYY,L);
            xm=mod(x3,L);ym=mod(y3,L);  
    
            urhs4 = interp2(xx,yy,uinterp,xm,ym,inter,0);
            vrhs4 = interp2(xx,yy,vinterp,xm,ym,inter,0); 
                %Instead of interping on a subset xm (a 256 grid of len L/3) of xx (a 768 grid of len L), 
                %we will interp on a finer set xm of xxcoarse (a 256 grid of len L):
                    %urhs4 = interp2(xxcoarse,yycoarse,uinterp,xm,ym,inter,0);
                    %vrhs4 = interp2(xxcoarse,yycoarse,vinterp,xm,ym,inter,0);
     %********************************************************************   
        
            dx = (deltat/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
            dy = (deltat/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
    
            pointsXXafter = pointsXX + dx;
            pointsYYafter = pointsYY + dy;
    
    %         figure;
    %         scatter(pointsXX,pointsYY,[],'filled')
    %         hold on;
    %         scatter(pointsXXafter,pointsYYafter)
    %         hold off;
    
    
    %         [G11, G12] = gradient(X,xss,yss);
    %         [G21, G22] = gradient(Y,xss,yss);
            G11 = (1/(2*h))*(pointsXXafter(2,3) - pointsXXafter(2,1));
            G12 = (1/(2*h))*(pointsXXafter(3,2) - pointsXXafter(1,2));
    
            G21 = (1/(2*h))*(pointsYYafter(2,3) - pointsYYafter(2,1));
            G22 = (1/(2*h))*(pointsYYafter(3,2) - pointsYYafter(1,2));
            %compute Cauchy-Green tensor
            CG11 = G11.*G11 + G21.*G21 ;
            CG12 = G11.*G12 + G21.*G22 ;
            
            CG21 = G12.*G11 + G22.*G21 ;
            CG22 = G12.*G12 + G22.*G22 ;
            
            Tr=CG11+CG22;
            Del=CG11.*CG22-CG12.*CG21;
            lam=Tr/2+sqrt(Tr.^2-4*Del)/2;
            %dle=log(lam)/(2)
            %dle=log(lam)/(2*INTTIME);            %Use this for net difference
            dle=log(lam)/(2*deltat);   %Use this for incremental diff


            %Now save the dle history
            if caseSelect == 1
                FTLEHistoryHighPoints(i) = dle;
            elseif caseSelect == 2
                FTLEHistoryMedPoints(i) = dle;
            else
                FTLEHistoryLowPoints(i) = dle;
            end
    
            
            %We want to compare this FTLE 4 particle to Strain 4 particle:

           

%              if caseSelect == 1
%                 rowIndex = rowIndexHigh;
%                 colIndex = colIndexHigh;
%             elseif caseSelect == 2
%                 rowIndex = rowIndexMed;
%                 colIndex = colIndexMed;
%             else
%                 rowIndex = rowIndexLow;
%                 colIndex = colIndexLow;
%             end

%             if caseSelect == 1
% 
%                 figure;
%                 quiver(u1Data,v1Data)
%     
%                rowIndexHigh
%                colIndexHigh
% 
%                u1Data(rowIndexHigh,colIndexHigh)
% 
%                v1Data(rowIndexHigh,colIndexHigh)
%             end

%             A11 = (256/2)*(u1Data(rowIndex,colIndex+1) - u1Data(rowIndex,colIndex-1));
%             A12 = (256/2)*(u1Data(rowIndex+1,colIndex) - u1Data(rowIndex-1,colIndex));
%     
%             A21 = (256/2)*(v1Data(rowIndex,colIndex+1) - v1Data(rowIndex,colIndex-1));
%             A22 = (256/2)*(v1Data(rowIndex+1,colIndex) - v1Data(rowIndex-1,colIndex));

%             A11 = (256/2)*(u1Data(rowIndex+1,colIndex) - u1Data(rowIndex-1,colIndex));
%             A12 = (256/2)*(u1Data(rowIndex,colIndex+1) - u1Data(rowIndex,colIndex-1));
%     
%             A21 = (256/2)*(v1Data(rowIndex+1,colIndex) - v1Data(rowIndex-1,colIndex));
%             A22 = (256/2)*(v1Data(rowIndex,colIndex+1) - v1Data(rowIndex,colIndex-1));
            
%             pointsXX
% 
%             pointsYY

            u1points = interp2(xa,ya,u1Data,pointsXX,pointsYY,inter,0);
            v1points = interp2(xa,ya,v1Data,pointsXX,pointsYY,inter,0);


            A11 = (1/(2*h))*(u1points(2,2+1) - u1points(2,2-1));
            A12 = (1/(2*h))*(u1points(2+1,2) - u1points(2-1,2));
    
            A21 = (1/(2*h))*(v1points(2,2+1) - v1points(2,2-1));
            A22 = (1/(2*h))*(v1points(2+1,2) - v1points(2-1,2));


            A = [A11, A12; A21, A22];

            ASym = .5*(A + A');

            lamA = eig(ASym);
            LAM = lamA(2);
            %strain=log(LAM); %Strain should be log(LAM) ?
            strain = LAM;     %Then it is e^strain that should = FTLE

%             mustBeReal(LAM)

            if caseSelect == 1
                StrainHistoryHigh(i) = strain;
            elseif caseSelect == 2
                StrainHistoryMed(i) = strain;
            else
                StrainHistoryLow(i) = strain;
            end

        end



        
        %*******************************************************************
        %*******************************************************************


        %Now I wish to focus on a particle trajectory (like Low FTLE), and
        %see how the nearby points evolve. Take a circle of nearby points
        %and plot the following ellipses as a movie:

        %First plot initial circle in red
        scatter(vxHistoryPlot(:,1),vyHistoryPlot(:,1),2,'filled','r') 
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        pbaspect([1 1 1])
        %Then plot the current ellipse.
%         hold on;
% 
% 
%         hold off;

        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    

    dle_2d;
    
    FTLEHigh = dle(rowIndexHigh,colIndexHigh);
    FTLEMed = dle(rowIndexMed,colIndexMed);
    FTLELow = dle(rowIndexLow,colIndexLow);


%     %Looks like strain needs to be shifted. 
%     for index3 = 1:no_of_steps-1
%         StrainHistoryHigh(index3) = StrainHistoryHigh(index3 + 1) ;
%         StrainHistoryMed(index3) = StrainHistoryMed(index3 + 1) ;
%         StrainHistoryLow(index3) = StrainHistoryLow(index3 + 1) ;
%     end

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
    plot(FTLEHistoryHigh,'k-')
    hold on;
    plot(FTLEHistoryHighPoints,'b--')
    hold on;
    plot(StrainHistoryHigh,'r:')
    hold on;
    plot(FTLEHigh*ones(1,no_of_steps),'-.')
    hold on;
    plot(FTLENetHistoryHigh,'LineWidth',2)
    hold off;
    legend('FTLE increment','FTLE 4 Points','Strain 4 Points','Net FTLE','Net FTLE History')
    title('FTLE High')

    subplot(3,1,2)
    plot(FTLEHistoryMed,'k-') 
    hold on;
    plot(FTLEHistoryMedPoints,'b--')
    hold on;
    plot(StrainHistoryMed,'r:')
    hold on;
    plot(FTLEMed*ones(1,no_of_steps),'-.')
    hold on;
    plot(FTLENetHistoryMed,'LineWidth',2)
    hold off;
    legend('FTLE increment','FTLE 4 Points','Strain 4 Points','Net FTLE','Net FTLE History')
    title('FTLE Med')

    subplot(3,1,3)
    plot(FTLEHistoryLow,'k-')
    hold on;
    plot(FTLEHistoryLowPoints,'b--')
    hold on;
    plot(StrainHistoryLow,'r:')
    hold on;
    plot(FTLELow*ones(1,no_of_steps),'-.')
    hold on;
    plot(FTLENetHistoryLow,'LineWidth',2)
    hold off;
    legend('FTLE increment','FTLE 4 Points','Strain 4 Points','Net FTLE','Net FTLE History')
    title('FTLE Low')

    

    



%Renew base time    
    start_time=start_time+BASESTEP;
end

%ss=strcat('save FTLEBTurb.mat FTLE');
%eval(ss);

