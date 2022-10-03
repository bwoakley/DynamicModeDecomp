function data_adv(XS,YS,xx,yy,varargin)
global domain resolution Nx Ny dxx
global s_dir rd_dir
global inter
global comp_type noise_type pu_type
global Period INTTIME CaseSTART N_Basetime dt pt nsteps nstep_outs ft n_cut
global periodic

nargin=numel(varargin);
if(nargin==2)
    xs=varargin{1};
    ys=varargin{2};
elseif(nargin==6)
    xz=varargin{1};
    yz=varargin{2};
    Mlen=varargin{3};
    file_num=varargin{4};
    flag=varargin{5};
    flagt=varargin{6};
else
;%nargin=0; this is for check filament
end

for ti=1:N_Basetime
    start_time=(ti-1)*pt+CaseSTART;
    ti
    xnow = XS;
    ynow = YS;
    index1=-10000;
    %if dir<0
    %    integ_time=min(start_time,INTTIME);
    %else
        integ_time=INTTIME;    
    %end
    no_of_steps = integ_time/abs(dt);
    
    for timeunit = flagt:no_of_steps
        i=timeunit;
        current_time = start_time + (i-1)*dt;   
        x0 = xnow;
        y0 = ynow;
        t0 = current_time;

        index1n = floor(t0/ft) + 1;
        if index1n~=index1
            index1 = index1n
            index2 = index1 + 1
            checkload;
        end    
        uinterp = u1 + (t0 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t0 - (index1-1)*ft)*(v2-v1)/ft;
        extraptolinear;

        urhs1 = interpn(xx,yy,UUU,x0,y0,inter,0);
        vrhs1 = interpn(xx,yy,VVV,x0,y0,inter,0);
        x1 = x0 + urhs1*dt/2;
        y1 = y0 + vrhs1*dt/2;
        t1 = t0 + dt/2;
        
        index1n = floor(t1/ft) + 1;
        if index1n~=index1
            index1 = index1n
            index2 = index1 + 1
            checkload;
        end 
        uinterp = u1 + (t1 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t1 - (index1-1)*ft)*(v2-v1)/ft;
        extraptolinear;
        
        urhs2 = interpn(xx,yy,UUU,x1,y1,inter,0);
        vrhs2 = interpn(xx,yy,VVV,x1,y1,inter,0);
        x2 = x0 + urhs2*dt/2;
        y2 = y0 + vrhs2*dt/2;
        t2 = t0 + dt/2;   
             
        urhs3 = interpn(xx,yy,UUU,x2,y2,inter,0);
        vrhs3 = interpn(xx,yy,VVV,x2,y2,inter,0);   
        x3 = x0 + urhs3*dt;
        y3 = y0 + vrhs3*dt;
        t3 = t0 + dt;
        
        index1n = floor(t3/ft) + 1;
        if index1n~=index1
            index1 = index1n
            index2 = index1 + 1
            checkload;
        end 
        uinterp = u1 + (t3 - (index1-1)*ft)*(u2-u1)/ft;
        vinterp = v1 + (t3 - (index1-1)*ft)*(v2-v1)/ft;
        extraptolinear;
        
        urhs4 = interpn(xx,yy,UUU,x3,y3,inter,0);
        vrhs4 = interpn(xx,yy,VVV,x3,y3,inter,0);   
        xnow = xnow+(dt/6)*(urhs1 + 2*urhs2 + 2*urhs3 + urhs4);
        ynow = ynow+(dt/6)*(vrhs1 + 2*vrhs2 + 2*vrhs3 + vrhs4);
            
        if(strcmp(comp_type,'fields'))
            if(timeunit==1),low_2d; end;
        else
            [W1,W2]=noisetype(numel(xnow),noise_type); 
            xnow=xnow+W1;
            ynow=ynow+W2;
            if(periodic(1))
                xnow=mod(xnow-domain(1,1),domain(1,2)-domain(1,1))+domain(1,1);
            end
            if(periodic(2))
                ynow=mod(ynow-domain(2,1),domain(2,2)-domain(2,1))+domain(2,1);
            else
                ynow(ynow>domain(2,2))=domain(2,2);
                ynow(ynow<domain(2,1))=domain(2,1);
            end
            if(mod(flag,nsteps)==0)
                if(strcmp(pu_type,'GPU'))
                    xtemp=gather(xnow(:)); ytemp=gather(ynow(:));
                else
                    xtemp=xnow(:); ytemp=ynow(:);
                end
                flag1=ceil(flag/nsteps)+1;
                
                xz(:,flag1)=xtemp(:);        % Tracer array at field times
                yz(:,flag1)=ytemp(:);
                if(mod(flag1,n_cut)==0)
                    xb=[]; yb=[];
                    xb{1}=xz(1:Mlen(1),:);
                    yb{1}=yz(1:Mlen(1),:);
                    a=Mlen(1);
                    if(size(Mlen,2)~=1)           %Separate the tracer sets
                        for ii=2:size(Mlen,2)
                            xb{ii}=xz(a+1:a+Mlen(ii),:);
                            yb{ii}=yz(a+1:a+Mlen(ii),:);
                            a=a+Mlen(ii);
                        end
                    end
                    savefile=strcat(s_dir,'final_positions',noise_type,num2str(file_num),'.mat');
                    sss=strcat('save',savefile,' xb yb Mlen');
                    eval(sss);
                    file_num=file_num+1;
                    xz=zeros(numel(XS),n_cut);
                    yz=zeros(numel(XS),n_cut);
                end
            end
            flag=flag+1;
        end
    end
    if(strcmp(comp_type,'fields'))
        meso_2d;
        dle_2d;
    end
end
end


function checkload
if(index1>9999)
    filename1 = strcat(rd_dir,'bin',num2str(index1));
elseif(index1>999)
    filename1 = strcat(rd_dir,'bin0',num2str(index1));
elseif(index1>99)
    filename1 = strcat(rd_dir,'bin00',num2str(index1));
elseif(index1>9)
    filename1 = strcat(rd_dir,'bin000',num2str(index1));
else
    filename1 = strcat(rd_dir,'bin0000',num2str(index1));
end

if(index2>9999)
    filename2 = strcat(rd_dir,'bin',num2str(index2));
elseif(index2>999)
    filename2 = strcat(rd_dir,'bin0',num2str(index2));
elseif(index2>99)
    filename2 = strcat(rd_dir,'bin00',num2str(index2));
elseif(index2>9)
    filename2 = strcat(rd_dir,'bin000',num2str(index2));
else
    filename2 = strcat(rd_dir,'bin0000',num2str(index2));
end

fid=fopen(filename1,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);

data_0=reshape(data,Nx,Ny,3);
u1=data_0(:,:,1);
v1=data_0(:,:,2);
w1=data_0(:,:,3);

fid=fopen(filename2,'rb');
data=fread(fid,[1 1],'*float');
data=fread(fid,[1 inf],'*double');
fclose(fid);

data_0=reshape(data,Nx,Ny,3);
u2=data_0(:,:,1);
v2=data_0(:,:,2);
end


function extraptolinear
for ii=1:3
    for jj=1:3
    UUU([1:Nx]+(ii-1)*Nx,[1:Ny]+(jj-1)*Ny,:)=uinterp;
    VVV([1:Nx]+(ii-1)*Nx,[1:Ny]+(jj-1)*Ny,:)=vinterp;
    end
end
end