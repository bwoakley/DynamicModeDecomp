Nx=144; Ny=64;
AI=zeros(Nx,Ny,20);
addpath ../
s_dir=strcat(' ../../data_revs/'); %save tracer directory
for t0=1:20
loader = strcat('load ',s_dir,'FP_set_ftcs_full_edit3_',num2str(t0),'.mat',' AA3');
saver = strcat('save ',s_dir,'FP_set_ftcs_full_1_',num2str(t0),'.mat',' AA3');

eval(loader)
AA=AA3;
ind1=find(AA==-100);
ind2=find(AA==-50);
ind3=find(AA==0);
ind4=find(AA==125);

BB=0.*AA;
BB(ind1)=1; 
BB(ind2)=-1; 
BB(ind3)=2;
BB(ind4)=0;

AA3=BB;

%AA=AA3;
figure(1);clf;
%subplot(2,1,1)
pcolor(reshape(AA3,Ny,Nx)); shading flat; drawnow
daspect([1 1 1])
Temp_geo=reshape(AA,Ny,Nx);
drawnow
pause
flag=0;
while flag==1
    prompt='change something?'
    [x_coord y_coord] = ginput();
    x_extra=[];
    if(isempty(x_coord))
        flag=0;
    else
    x_extra(:,1)=floor(x_coord);
    x_extra(:,2)=floor(y_coord);
    I=Temp_geo;
    prompt='what value?'
    w = waitforbuttonpress;
if w == 0
    disp('Button click')
    [x_coord2 y_coord2] = ginput(1);
   x_2=floor(x_coord2);
   y_2=floor(y_coord2);
   n_val=I(y_2,x_2)
else
    disp('Key press')
    n_val=input(prompt);
end


    end
    if(~isempty(x_extra))
        x2=x_extra;
        for jj=1:size(x_extra,1)
            I = flood2(I,n_val,x_extra(jj,2),x_extra(jj,1));
        end
    end
    Temp_geo=I;
 
    title(num2str(t0))    
    
    subplot(2,1,2)
pcolor(Temp_geo); shading flat; drawnow
daspect([1 1 1])
colorbar
    
end

% AA=0.*Temp_geo;
% AA(Temp_geo==n_val)=1;
% subplot(2,1,2)
% pcolor(AA); shading flat; drawnow
% daspect([1 1 1])
% colorbar

AA3=Temp_geo;

FTCS(:,:,t0)=AA3';
eval(saver)

end


save ftcs_sets.mat FTCS

