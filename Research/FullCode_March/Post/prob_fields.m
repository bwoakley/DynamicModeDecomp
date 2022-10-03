%% This post-processing subroutine is designed to work directly with the outputs of the main_code;
r_dir=strcat(' ../../data_revs/'); %read tracer directory
s_dir=strcat(' ../../data_revs/'); %save tracer directory

loader=strcat('load',r_dir,'globs.mat');
eval(loader);

prompt='Computing (FTE) or (FP)? ';
field_type=input(prompt);

prompt='Base time: ';
ti=input(prompt);
prompt='End time: ';
tf=input(prompt);

prompt='How many periods? ';
nP=input(prompt);

%%
periodic=[1 0]; %Should be editted
Nx=Nx/4; Ny=Ny/4; %This part just needs to be editted to fix memory
filelast1=0;
filelast2=0;
for t0=ti:tf
    save_fp = strcat('save ',s_dir,'FP_set_',num2str(t0),'.mat Ab Di1 Br B0 Nx Ny');
    save_fte = strcat('save ',s_dir,'FTE_field',num2str(t0),'.mat FTE x0 y0 xr yr');

    t0
    file_ind=ceil(t0/n_cut);
    flag=mod(t0,n_cut); 
    if(flag==0);
        flag=n_cut;
    end
    load1=strcat('load ',r_dir,'final_positions',noise_type,num2str(file_ind),'.mat');%[xb,yb,Mlen]
    
    

    if(filelast1~=file_ind) 
        eval(load1);        
        xz0(1:Mlen(1),:)=xb{1};               %Copy xz back over
        yz0(1:Mlen(1),:)=yb{1};
        a=Mlen(1);
        if(size(Mlen,2)~=1)           %Separate the tracer sets
            for ii=2:size(Mlen,2)
                xz0(a+1:a+Mlen(ii),:)=xb{ii};
                yz0(a+1:a+Mlen(ii),:)=yb{ii};
                a=a+Mlen(ii);
            end
        end
        filelast1=file_ind;
    end
    xnow=xz0(:,flag);
    ynow=yz0(:,flag);
    
    if(periodic(1))
        xnow=mod(xnow-domain(1,1),diff(domain(1,:)))+domain(1,1);
    end
    if(periodic(2))
        ynow=mod(ynow-domain(2,1),diff(domain(2,:)))+domain(2,1);
    else
        ynow(ynow>domain(2,2))=domain(2,2);
        ynow(ynow<domain(2,1))=domain(2,1);
    end

    x0=1+round((((Nx-1)/diff(domain(1,:)))*(xnow-domain(1,1))));
    y0=1+round((((Ny-1)/diff(domain(2,:)))*(ynow-domain(2,1))));
    B0=Nx*(y0-1)+x0;

    
    
    file_ind=ceil((t0+(nP*20))/n_cut);
    flag=mod(t0+ceil(nP*20),n_cut);
    if(flag==0);
        flag=n_cut;
    end
    load2=strcat('load ',r_dir,'final_positions',noise_type,num2str(file_ind),'.mat');%[xb,yb,Mlen]
    
    if(filelast2~=file_ind)
        eval(load2);                
        xzf(1:Mlen(1),:)=xb{1};               %Copy xz back over
        yzf(1:Mlen(1),:)=yb{1};
        a=Mlen(1);
        if(size(Mlen,2)~=1)           %Separate the tracer sets
            for ii=2:size(Mlen,2)
                xzf(a+1:a+Mlen(ii),:)=xb{ii};
                yzf(a+1:a+Mlen(ii),:)=yb{ii};
                a=a+Mlen(ii);
            end
        end
        filelast2=file_ind;
    end
    xnow=xzf(:,flag);
    ynow=yzf(:,flag);

    if(periodic(1))
        xnow=mod(xnow-domain(1,1),diff(domain(1,:)))+domain(1,1);
    end
    if(periodic(2))
        ynow=mod(ynow-domain(2,1),diff(domain(2,:)))+domain(2,1);
    else
        ynow(ynow>domain(2,2))=domain(2,2);
        ynow(ynow<domain(2,1))=domain(2,1);
    end

    xr=1+round((((Nx-1)/diff(domain(1,:)))*(xnow-domain(1,1))));
    yr=1+round((((Ny-1)/diff(domain(2,:)))*(ynow-domain(2,1))));
    Br=Nx*(yr-1)+xr;

    
    
switch field_type
    case 'FTE'
        FTEm=zeros(1,Nx*Ny);  
        for i=1:Nx*Ny
            ind=find(B0==i);
            delta=zeros(1,Nx*Ny);
            for j=1:length(ind)
                delta(Br(ind(j)))=delta(Br(ind(j)))+1;
            end
            delta=delta/length(ind);
            I=find(delta);
            FTEm(i)=-1*dot(delta(I),log(delta(I)));
        end
        FTE=reshape(FTEm,Nx,Ny);
        eval(save_fte)
        
    case 'FP'
        B=1:Nx*Ny;
        Di1=zeros(length(B(:)),1);
        Ab=zeros(length(B(:)));
        for i=1:length(B)
            ind1=find(B0==B(i));
            if(~isempty(ind1))
                Di1(B(i))=Di1(B(i))+numel(ind1);
                Bj=unique(sort(Br(ind1)));
                for j=1:numel(Bj)
                    ind2=find(Br(ind1)==Bj(j));
                    if(~isempty(ind2))
                        Ab(B(i),Bj(j))=Ab(B(i),Bj(j))+numel(ind2);
                    end    
                end
                Ab(B(i),:)=Ab(B(i),:)/numel(ind1);
            end
        end
        eval(save_fp)
    end
end