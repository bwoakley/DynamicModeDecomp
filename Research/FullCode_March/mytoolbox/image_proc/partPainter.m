Nx=128; Ny=128;
%FTCS=zeros(Nx,Ny,20);
addpath ../Desktop/Results_PV/Partitions/Turb/real_GEO
%s_dir=strcat(' ../../data_revs/'); %save tracer directory
s_dir=strcat(' ');
for t0=91:-1:1
loader1 = strcat('load ',s_dir,' time_eidt6_',num2str(t0),'.mat',' AA3');
%loader2 = strcat('load ',s_dir,' time_eidt6_',num2str(t0+1),'.mat',' AA3');
saver = strcat('save ',s_dir,' time_eidt6_',num2str(t0),'.mat',' AA3');

%eval(loader2)
%BB=AA3;
eval(loader1)

AA=AA3;
f1=figure(1);clf;
colormap('colorcube')
s1=subplot(2,2,1);
pcolor(reshape(AA,Ny,Nx)); shading flat; drawnow
daspect([1 1 1])
 s2=subplot(2,2,2);
 pcolor(reshape(AA,Ny,Nx)); shading flat; drawnow
daspect([1 1 1])
s3=subplot(2,2,3);
pcolor(reshape(AA,Ny,Nx)); shading flat; drawnow
daspect([1 1 1])
s4=subplot(2,2,4);
pcolor(reshape(AA,Ny,Nx)); shading flat; drawnow
daspect([1 1 1])


xstart=.06;
xwidth=.4;
ywidth=.5;

set(f1,'position',[10 10 800 800]);
set(s1,'position',[xstart+0*(xwidth+.0) .46 xwidth ywidth])
set(s2,'position',[xstart+1*(xwidth+.0) .46 xwidth ywidth])
set(s3,'position',[xstart+0*(xwidth+.05) .06 xwidth ywidth])
set(s4,'position',[xstart+1*(xwidth+.0) .06 xwidth ywidth])


% subplot(2,1,2)
% pcolor(field); shading interp; drawnow
% daspect([1 1 1])
Temp_geo=reshape(AA,Ny,Nx);
drawnow
flag=1;
%pause()

% prompt='what value?'
% n_val=input(prompt);
I=Temp_geo;
while flag==1
    
set(f1,'position',[10 10 800 800]);
set(s1,'position',[xstart+0*(xwidth+.01) .46 xwidth ywidth])
set(s2,'position',[xstart+1*(xwidth+.0) .46 xwidth ywidth])
set(s3,'position',[xstart+0*(xwidth+.01) .06 xwidth ywidth])
set(s4,'position',[xstart+1*(xwidth+.0) .06 xwidth ywidth])
    
    prompt='change something?'
    
    
    w = waitforbuttonpress;
if w ~= 0
    disp('Button click')
%     pcolor(Temp_geo); shading flat; drawnow
%     daspect([1 1 1])
    
    
    s1=subplot(2,2,1);
 pcolor(Temp_geo); shading flat; drawnow
daspect([1 1 1])
s2=subplot(2,2,2);
pcolor(Temp_geo); shading flat; drawnow
daspect([1 1 1])
s3=subplot(2,2,3);
pcolor(Temp_geo); shading flat; drawnow
daspect([1 1 1])
s4=subplot(2,2,4);
pcolor(Temp_geo); shading flat; drawnow
daspect([1 1 1])
    

else
    Temp_geo=I;
    [x_coord y_coord] = ginput(1);
    x_extra=[];
    if(isempty(x_coord))
        flag=0;
    else
    x_extra(:,1)=floor(x_coord);
    x_extra(:,2)=floor(y_coord);
    I=Temp_geo;
    
n_val=40;
    end
    if(~isempty(x_extra))
        x2=x_extra;
        for jj=1:size(x_extra,1)
             I = flood2(I,n_val,x_extra(jj,2),x_extra(jj,1));
%        ind=find(Temp_geo==Temp_geo(x_extra(jj,2),x_extra(jj,1)));
%        I(ind)=n_val;
        end
    end
 
    title('good?')    
    
% pcolor(I); shading flat; drawnow
% daspect([1 1 1])

s1=subplot(2,2,1);
 pcolor(I); shading flat; drawnow
daspect([1 1 1])
s2=subplot(2,2,2);
pcolor(I); shading flat; drawnow
daspect([1 1 1])
s3=subplot(2,2,3);
pcolor(I); shading flat; drawnow
daspect([1 1 1])
s4=subplot(2,2,4);
pcolor(I); shading flat; drawnow
daspect([1 1 1])


end

AA3=Temp_geo;
eval(saver)
end



end
