%This code uses the extrema2 subroutine

% X=gather(xnow(:,1:Nx));
% Y=gather(ynow(:,1:Nx));

X=gather(xnow(1:Nx*Ny)); Y=gather(ynow(1:Nx*Ny));
X=reshape(X,Ny,Nx); Y=reshape(Y,Ny,Nx);



xs=gather(xs);
ys=gather(ys);

[ux uy]=gradient(X,xs,ys);
[vx vy]=gradient(Y,xs,ys);

Q=((ux-vy).^2+(vx+uy).^2-(vx-uy).^2);
Q0=std(Q(:));
Qnow=Q./Q0;

%The outputs need to be the OW field, and the elliptic centers
%(PoincareSections) to be used for the geodesic elliptic extraction.
[xymax,smax,xymin,smin]=extrema2(Qnow');
min_val=min(Qnow(:));
ind=find(xymin<=-min_val);
smin=smin(ind);

y_loc=floor(smin./(Nx));
ind=find(y_loc==0);
y_loc=y_loc+1;
x_loc=smin-Nx*(y_loc-1);
x_loc=x_loc+1;
yyloc=ys(y_loc);
xxloc=xs(x_loc);

%The code now accurately obtains the highly ellptic regions
%Should figure out how to automate this part...
poincareSection = struct('endPosition',{},'numPoints',{},'orbitMaxLength',{});
for i=1:size(smin,1)
    if(yyloc(i)>0)
        poincareSection(i).endPosition = [xxloc(i),yyloc(i)+.75;xxloc(i)-1.8,yyloc(i)+.75];
    else
        poincareSection(i).endPosition = [xxloc(i),yyloc(i)-.75;xxloc(i)+1.8,yyloc(i)-.75];
    end
end

M=size(smin,1);
savefile=strcat(s_dir,'ow_field',num2str(ti),'.mat');
sss=strcat(' save',savefile,' Qnow poincareSection M');
eval(sss);

% subplot(1,3,1)
% pcolor(xs,ys,Qnow); shading interp; shg
% hold on
% plot(xxloc,yyloc,'kp'); shg
