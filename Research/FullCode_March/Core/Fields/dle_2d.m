function dle_2d
%This code requires the eig_array subroutine
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

 gradF11 = G11;
 gradF12 = G12;
 gradF21 = G21;
 gradF22 = G22;
    
 gradF11 = reshape(gradF11,prod(double(resolution)),1);
 gradF12 = reshape(gradF12,prod(double(resolution)),1);
 gradF21 = reshape(gradF21,prod(double(resolution)),1);
 gradF22 = reshape(gradF22,prod(double(resolution)),1);

cgStrainMainGrid(:,1) = gradF11.^2 + gradF21.^2;
cgStrainMainGrid(:,2) = gradF11.*gradF12 + gradF21.*gradF22;
cgStrainMainGrid(:,3) = gradF12.^2 + gradF22.^2;

nRows = size(cgStrainMainGrid,1);
[cgStrainV,cgStrainD] = arrayfun(@(x11,x12,x22)eig_array(x11,x12,x22),cgStrainMainGrid(:,1),cgStrainMainGrid(:,2),cgStrainMainGrid(:,3),'uniformOutput',false);
cgStrainV = cell2mat(cgStrainV);
cgStrainD = cell2mat(cgStrainD);
            
nRows = size(cgStrainMainGrid,1)
cgStrain = arrayfun(@(idx)[cgStrainMainGrid(idx,1),cgStrainMainGrid(idx,2);cgStrainMainGrid(idx,2),cgStrainMainGrid(idx,3)],1:nRows,'uniformOutput',false);
cgStrain = cell2mat(cgStrain);
cgStrain = reshape(cgStrain,[2 2 nRows]);

idx=cgStrainD(:,2)<1;
n=sum(idx);
cgStrainD(~idx,1)=1./cgStrainD(~idx,2);
cgStrainD(idx,:)=nan;
cgStrainV(idx,:)=nan;

negIdx=any(cgStrainD<=0,2);
cgStrainD(negIdx,:)=nan;
cgStrainV(negIdx,:)=nan;
cgStrain(:,:,negIdx)=nan;

cg2=reshape(cgStrainD(:,2),fliplr(resolution));
dle=.5*log(cg2)./INTTIME;

% subplot(1,3,2)
% pcolor(xs,ys,dle); shading interp

savefile=strcat(s_dir,'ftle_field',num2str(ti),'.mat');
sss=strcat('save',savefile,' cgStrainV cgStrainD cgStrain dle domain resolution');
eval(sss)
%By outputting the eigenvector/value-fields the elliptic geodesic
%structures can be more quickly computed.
end

function [v,d] = eig_array(x11,x12,x22)

[v,d] = eig_custom([x11,x12;x12,x22]);
d = transpose(diag(d));
v = reshape(v,1,4);

function [v,d] = eig_custom(a)

d(2,2) = .5*trace(a) + sqrt(.25*trace(a)^2 - det(a));
d(1,1) = .5*trace(a) - sqrt(.25*trace(a)^2 - det(a));
if any(imag(d(:)))
    warning([mfilename,':complexEigenvalue'],['Complex eigenvalue: ',num2str(d([1,4]))])
end

v(1,2) = -a(1,2)/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);
v(2,2) = (a(1,1) - d(2,2))/sqrt(a(1,2)^2 + (a(1,1) - d(2,2))^2);

v(1,1) = v(2,2);
v(2,1) = -v(1,2);