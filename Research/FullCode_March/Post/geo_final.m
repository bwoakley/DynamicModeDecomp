clear all; clc; close all;
k=1;
r_dir=strcat(' fields/'); %read directory
s_dir=strcat(' fields/'); %save directory

for ti=1
    loadfile=strcat(r_dir,'geo_poss_chk',num2str(ti),'.mat');
    sss=strcat('load',loadfile);
    eval(sss);
   resolution=[144*4,64*4];
   domain=[[-2.5*pi 2.5*pi];[-3 3]];
    Nx=resolution(1);
    Ny=resolution(2);
    elx=linspace(domain(1,1),domain(1,2),k*Nx+1); elx=elx(1:end-1);
    ely=linspace(domain(2,1),domain(2,2),k*Ny+1); ely=ely(1:end-1);
    [yelx xelx]=ndgrid(ely,elx);
    xs=linspace(domain(1,1)-diff(domain(1,:)),domain(1,2)+diff(domain(1,:)),3*Nx+1);xs = xs(1:end-1);
    [yst, xst] = ndgrid(ely,xs);
    
    for jj=1:MM
        geodesX{jj+MM}(:,1)=geodesX{jj}(:,1)+5*pi;
        geodesX{jj+MM}(:,2)=geodesX{jj}(:,2)+0;
        geodesY{jj+MM}(:,1)=geodesY{jj}(:,1)+5*pi;
        geodesY{jj+MM}(:,2)=geodesY{jj}(:,2)+0;
    end
    for jj=1:MM
        geodesX{jj+2*MM}(:,1)=geodesX{jj}(:,1)-5*pi;
        geodesX{jj+2*MM}(:,2)=geodesX{jj}(:,2)+0;
        geodesY{jj+2*MM}(:,1)=geodesY{jj}(:,1)-5*pi;
        geodesY{jj+2*MM}(:,2)=geodesY{jj}(:,2)+0;
    end    
   
    ellip_in=[];
    Temp=0.*xst; 
    %cfinal=0.*xnow(:); 
    for jj=1:3*MM
        in=inpolygon(xst,yst,geodesX{jj}(:,1),geodesX{jj}(:,2));
        %in=inpolygon(xnow(:),ynow(:),geodesX{jj}(:,1),geodesX{jj}(:,2));
        ind=find(in);
        ellip_in=[ellip_in; ind];
        
        in=inpolygon(xst,yst,geodesY{jj}(:,1),geodesY{jj}(:,2));
        %in=inpolygon(xnow(:),ynow(:),geodesX{jj}(:,1),geodesX{jj}(:,2));
        ind=find(in);
        ellip_in=[ellip_in; ind];
    end
    
    %cfinal(ellip_in)=1; %Indicate that they are in the elliptic region for that timestep 
    Temp(ellip_in)=1;
    B1=0.*xelx;
    A1=Temp(:,1:Nx);
    A2=Temp(:,(Nx+1):2*Nx);
    A3=Temp(:,(2*Nx+1):end);
    B1=A1+A2+A3;
    ind1=find(B1>0);
    B1(ind1)=1;
    
    figure(1); clf;
    plot(xelx(ind1),yelx(ind1),'.k'); hold on
    hClosedOrbitsEtaPos = arrayfun(@(i)plot(geodesX{i}(:,1),geodesX{i}(:,2)),1:3*MM);
    hClosedOrbitsEtaPos = arrayfun(@(i)plot(geodesY{i}(:,1),geodesY{i}(:,2)),1:3*MM);
    drawnow
    pause(.1)
    
    MM=3*MM;
    savefile=strcat(s_dir,'geo_efield_',num2str(ti),'.mat');
    sss=strcat('save',savefile,' B1');
    eval(sss);
end