function meso_2d
%This code requires the hypergraph plotting subroutine
X=gather(xnow); Y=gather(ynow);

xa=reshape(X(1:Nx*Ny),Ny,Nx); ya=reshape(Y(1:Nx*Ny),Ny,Nx);
xb=reshape(X(Nx*Ny+1:2*Nx*Ny),Ny,Nx); yb=reshape(Y(Nx*Ny+1:2*Nx*Ny),Ny,Nx);
xc=reshape(X(2*Nx*Ny+1:3*Nx*Ny),Ny,Nx); yc=reshape(Y(2*Nx*Ny+1:3*Nx*Ny),Ny,Nx);
xd=reshape(X(3*Nx*Ny+1:end),Ny,Nx); yd=reshape(Y(3*Nx*Ny+1:end),Ny,Nx);

% xa=X(:,1:Nx); ya=Y(:,1:Nx);
% xb=X(:,Nx+1:2*Nx); yb=Y(:,Nx+1:2*Nx);
% xc=X(:,2*Nx+1:3*Nx); yc=Y(:,2*Nx+1:3*Nx);
% xd=X(:,3*Nx+1:end); yd=Y(:,3*Nx+1:end);
    
G11=(xa-xb)/(2*dxx);
G12=(xc-xd)/(2*dxx);
G21=(ya-yb)/(2*dxx);
G22=(yc-yd)/(2*dxx);
    
ux=(G11-1)/INTTIME;
uy=G12/INTTIME;
vx=G21/INTTIME;
vy=(G22-1)/INTTIME;

DetVbar=ux.*vy-uy.*vx;
saver=strcat(s_dir,'meso_field',num2str(ti),'.mat');
sss=strcat(' save',saver,' DetVbar INTTIME');
eval(sss);

%subplot(1,3,3)
%fct_plot_hypergraph(xs,ys,DetVbar,INTTIME/Period,300)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author:   Sophie Loire
%           (sloire@engineering.ucsb.edu)
%           UCSB, Santa Barbara, CA 93106
%           for Hypergraph Project
%
% Date created 05/02/2014
% Last modified 11/24/2014
%
%   PURPOSE: Plot the MESO properties of a flow with an appropriate
%   red(Meso-Helical), green(Meso-elliptic) and blue(Meso-Hyperbolic)  colormap                
%   
%   INPUTS:
%       - z1, z2, coordinate for D
%       - D is the matrix containing the determinant of the gradient
%        of the average Lagrangian velocity
%       - T is the time scale used for the average Lagrangian velocity
%       - nbcolor is the number of colors in the colormap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function fct_plot_hypergraph(z1,z2,D,T, nbcolor)


%% Hypergraph colormap
cmap = zeros(nbcolor,3);
nb3 = floor(nbcolor/3);
ind = 0;
for i=1:nb3
    ind =ind +1;
    cmap(ind,1)=(i-1)/(nb3-1);
    cmap(ind,2)=(i-1)/(nb3-1);
    cmap(ind,3)=1;
end
for i=1:nb3
    ind =ind +1;
    cmap(ind,1)=(cos(pi*(i-1)/(nb3-1)))^2;
    cmap(ind,2)=1;
    cmap(ind,3)=(cos(pi*(i-1)/(nb3-1)))^2;
end
for i=1:nb3
    ind =ind +1;
    cmap(ind,1)=1;
    cmap(ind,2)=(nb3-i)/(nb3-1);
    cmap(ind,3)=(nb3-i)/(nb3-1);
end
colormap(cmap);
colorbar;

%% End colorbar

Ts= 4./T^2;
AA=(D>T);          % Meso Helical
BB=(D<0.);          % Meso Hyperbolic
CC=D;

%% Some rescaling and cutoffs
CC(AA)= Ts + log(1-Ts+D(AA));
CC(BB)= -log(1-D(BB));
CC(CC>2.*Ts)=2.*Ts;
CC(CC<-Ts)=-Ts;
CC(isnan(CC))=-(Ts+0.01);

%% Plotting

hh=pcolor(z1,z2,CC);
set(hh,'edgecolor','none');
caxis([-(Ts+0.01) 2*Ts]);
colormap(cmap);
colorbar;
end

