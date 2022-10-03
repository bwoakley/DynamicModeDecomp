clear all; clc; close all;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-----------Point to subroutines that will be needed-----------------------

global dxx Nx Ny domain resolution
global s_dir rd_dir
global inter
global comp_type noise_type pu_type
global Period INTTIME CaseSTART N_Basetime dt pt nsteps n_cut ft
global periodic
tic
% addpath ../mytoolbox/extrema/
% addpath ../mytoolbox/stbl/
% addpath ../mytoolbox/geo_tools/
% addpath ../mytoolbox/
inter = '*linear';
r_dir=strcat(' ../../data_revs/'); %read tracer directory
s_dir=strcat(' ../../data_revs/'); %save tracer directory
rd_dir=strcat('/media/philw/Titan/2d_turb/advcore/data/res/'); % read in the velocity fields data
sf_dir=strcat(' ../../data_revs/'); %save field directory
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------------------User input---------------------------------------------
prompt='Compute on what? ';
pu_type=input(prompt);
if(strcmp(pu_type,'GPU'))
    gpuDevice();
end

prompt='Computing (fields) or (tracers)? ';
comp_type=input(prompt);
      
prompt='What flow field? [Bick,Stand,DG,Turb] ';
field_case=input(prompt);

prompt='(eq) or (data)? ';
eq_type=input(prompt);   
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-------This section sets up the domain and file locations-----------------

switch field_case
    case 'Turb'
        Lx=2*pi; Ly=2*pi;   %Set domain length
        domain=[0 Lx;0 Ly];
        Nx=256; Ny=256;     %Field resolutions
        Period=1;           %If no period set as 1
    case 'Bick'
        Lx=5*pi; Ly=6;
        domain=[[-1 1]*Lx/2;[-1 1]*Ly/2];
        Nx=144*1; Ny=64*1;
        param.beta=.6144;
        param.e1=0.002;
        param.e2=0.3;
        param.Delta=sqrt(1-param.beta*3/2);
        param.c1=1/3*(1+param.Delta);
        param.c2=1/3*(1-param.Delta);
        param.k1=sqrt(6*param.c1);
        param.k2=sqrt(6*param.c2); 
        param.Om=param.k1*(param.c1-param.c2);
        Period=2*pi/param.Om;
        addpath ../flows/Bic/       
     case 'DG'
         Lx=2; Ly=2;
         domain=[0 Lx; 0 Ly];
         Nx=256; Ny=256;
         Period=1;
         addpath ../flows/DG/
    case 'Thin'
         Lx=4; Ly=4;
         domain=[[0 1]*Lx;[0 1]*Ly];
         Nx=256; Ny=256;
         Period=1;
         addpath ../flows/Thin/
     case 'Other'
%         Lx=; Ly=;
%         domain=[];
%         Nx=; Ny=;
%         Period=;
         addpath ../flows/other/
end
resolution=[Nx,Ny];
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%-----------------------Set other parameters-------------------------------

dir=1;                             %Direction of time
ft=.002*Period;                    %File separation (if using data)
%ft=1;
pt=0.05*Period;                    %Spacing between partitions
%pt=100;
INTTIME=11*Period+pt;                 %Integration time
nsteps=pt/ft;                      %How often to store the tracer locations
n_cut=25;
dt = dir*ft;
dxx=1e-6;
e_vec=linspace(.014,.016,10);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%----------------------------Setup the grids to use------------------------

switch comp_type
    case 'tracers'
        nstep_outs=n_cut*nsteps;    %How often to save the tracer histories
        prompt='Periodic?';
        periodic=input(prompt);   
        prompt='What type of noise?[N,G,L]';
        noise_type=input(prompt);
        prompt='Load from file? ';
        read_type=input(prompt);
        
        if(read_type)
            prompt='(initial) or (other)?';
            init_type=input(prompt);
            
            switch init_type
                case 'initial'
                    file_num=1; flagt=1; time_ind=0; flag=1;
                    loadfile=strcat(r_dir,'initial_positions.mat');%[XS,YS,Mlen]
                    sss=strcat('load',loadfile);
                    eval(sss);
                    
                    xz=zeros(numel(XS),n_cut);
                    yz=zeros(numel(XS),n_cut);
                    xz(:,1)=XS; yz(:,1)=YS;
                case 'other'
                    prompt='From what file index?';
                    time_ind=input(prompt);
                    flagt=(time_ind-1)*nsteps+1;
                    file_ind=ceil(time_ind/n_cut);
                    flag=mod(time_ind,n_cut); 
                    if(flag==0)
                        file_num=file_ind+1;
                        flag=n_cut;
                    else
                        file_num=file_ind;
                    end

                    loadfile=strcat(r_dir,'final_positions',noise_type,num2str(file_ind),'.mat');%[xb,yb,Mlen]
                    sss=strcat('load',loadfile);
                    eval(sss);
                    
                    xz(1:Mlen(1),:)=xb{1};               %Copy xz back over
                    yz(1:Mlen(1),:)=yb{1};
                    a=Mlen(1);
                    if(size(Mlen,2)~=1)           %Separate the tracer sets
                        for ii=2:size(Mlen,2)
                            xz(a+1:a+Mlen(ii),:)=xb{ii};
                            yz(a+1:a+Mlen(ii),:)=yb{ii};
                            a=a+Mlen(ii);
                        end
                    end 
            end
            if(strcmp(pu_type,'GPU'))
                XS=gpuArray(xz(:,flag)); % This should put only the correct time info, and not all
                YS=gpuArray(yz(:,flag)); % While still having read-in the entire xz history 
            else
                XS=xz(:,flag);
                YS=yz(:,flag);
            end
        else
            file_num=1; flag=1; time_ind=0; flagt=1;
            prompt='The square of tracers per grid:';
            k=input(prompt);
            xs=linspace(domain(1,1),domain(1,2),k*Nx+1);dx=xs(2)-xs(1); xs=xs(1:end-1)+dx/2;
            ys=linspace(domain(2,1),domain(2,2),k*Ny+1);dy=ys(2)-ys(1); ys=ys(1:end-1)+dy/2;
            [YS,XS]=ndgrid(ys,xs);
            xz=zeros(numel(XS),n_cut);
            yz=zeros(numel(XS),n_cut);
            xtemp=XS(:); ytemp=YS(:);
            xz(:,1)=xtemp(:); yz(:,1)=ytemp(:);
            Mlen=size(XS(:),1);
            if(strcmp(pu_type,'GPU'))
                XS=gpuArray(XS(:));
                YS=gpuArray(YS(:));
            else
                XS=XS(:);
                YS=YS(:);
            end  
        end
        CaseSTART=(time_ind-1)*(pt/Period);%Start at t=0, find LCS at different base times
        CaseEND=CaseSTART;        %End at t=1 period. End=Start for integration

    case 'fields'
        addpath Fields/
        prompt='Start at partition number?';
        file_num=input(prompt);
        prompt='End at partition number?';
        file_num2=input(prompt);
        
        xs=linspace(domain(1,1),domain(1,2),Nx+1)+dxx; xs=xs(1:end-1);
        ys=linspace(domain(2,1),domain(2,2),Ny+1); ys=ys(1:end-1);
        [YS1,XS1]=ndgrid(ys,xs);
        xs=linspace(domain(1,1),domain(1,2),Nx+1)-dxx; xs=xs(1:end-1);
        ys=linspace(domain(2,1),domain(2,2),Ny+1); ys=ys(1:end-1);
        [YS2,XS2]=ndgrid(ys,xs);
        xs=linspace(domain(1,1),domain(1,2),Nx+1); xs=xs(1:end-1);
        ys=linspace(domain(2,1),domain(2,2),Ny+1)+dxx; ys=ys(1:end-1);
        [YS3,XS3]=ndgrid(ys,xs);
        xs=linspace(domain(1,1),domain(1,2),Nx+1); xs=xs(1:end-1);
        ys=linspace(domain(2,1),domain(2,2),Ny+1)-dxx; ys=ys(1:end-1);
        [YS4,XS4]=ndgrid(ys,xs);
        YS=[YS1 YS2 YS3 YS4];
        XS=[XS1 XS2 XS3 XS4];
        if(strcmp(pu_type,'GPU'))
            XS=gpuArray(XS(:));
            YS=gpuArray(YS(:));
        else
            XS=XS(:);
            YS=YS(:);
        end
        noise_type='N';
        CaseSTART=(file_num-1)*(pt/Period);%Start at t=0, find LCS at different base times
        CaseEND=(file_num2-1)*(pt/Period); %End at t=1 period.
end
xx=linspace(domain(1,1)-Lx,domain(1,2)+Lx,Nx*3+1);xx=xx(1:end-1);
yy=linspace(domain(2,1)-Ly,domain(2,2)+Ly,Ny*3+1);yy=yy(1:end-1);
Diff=CaseEND-CaseSTART;          %Gap to compute fields  
N_Basetime=Diff/(pt/Period)+1;   %Number of partitions to compute
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%------------Start computing-----------------------------------------------


switch eq_type
    case 'data'
        if(strcmp(comp_type,'fields'))
            data_adv(XS,YS,xx,yy,xs,ys); %Use data to compute the fields
        else
            data_adv(XS,YS,xx,yy,xz,yz,Mlen,file_num,flag,flagt); %Use data to compute traj.
        end
    case 'eq'
        if(strcmp(comp_type,'fields'))
            eq_adv(XS,YS,xs,ys);  %Use equations to compute the fields
        else
            eq_adv(XS,YS,xz,yz,Mlen,file_num,flag,flagt); %Use equations to compute traj.
        end
end

%xs and YS make the computational domain, and I should split this
%computational task amoungst the cores of the GPU

if(strcmp(eq_type,'eq'))
    saver=strcat('save',s_dir,'globs.mat Nx Ny domain resolution dxx'...
    ,' s_dir comp_type noise_type dt pt n_cut Period INTTIME eq_type');
else
    saver=strcat('save',s_dir,'globs.mat Nx Ny domain resolution dxx'...
    ,'s_dir comp_type noise_type dt pt n_cut Period INTTIME eq_type rd_dir xx yy');
end
eval(saver);
% This is where I can append the post-processing subroutines.
toc

