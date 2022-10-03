% This is the main post-processing utility for tracers
global nP Nx Ny domain noise_type pt
global n_cut
tic
prompt='Compute on what? ';
pu_type=input(prompt);

prompt='(moms) or (resi)? ';
stat_type=input(prompt);

prompt='(N) (G) (L)? ';
noise_type=input(prompt);

rd_dir=' /media/philw/Thanos/fullcode_data/dg_data_test/';
r_dir=' /media/philw/Thanos/fullcode_data/dg_data_test/';%Fill this
s_dir=' /media/philw/Thanos/fullcode_data/dg_data_test/';

loader=strcat('load ',r_dir,'globs.mat');
eval(loader);

Narteng=INTTIME/(pt*n_cut);
Narteng=8;
Plang=n_cut;
Beltin=20;
nsteps=25;                      %How often to store the tracer locations

Lx=diff(domain(1,:));
Ly=diff(domain(2,:));
K=8;
xs=linspace(domain(1,1),domain(1,2),Nx+1); xs=xs(1:end-1);
ys=linspace(domain(2,1),domain(2,2),Ny+1); ys=ys(1:end-1);


switch stat_type
    case 'resi'

    prompt='(PS)(OW)(FTLE)(MESO)(GEO)(FTE)(AI)(FTCS)? ';
    res_type=input(prompt);

    for j=1:Narteng 
        xl=[];
        load_data=strcat('load ',rd_dir,'final_positions',noise_type,num2str(ti),'.mat');
        
        load_ps=strcat('load ',r_dir,'psfields_',num2str(ti+(j-1)*Plang),'.mat');
        load_ow=strcat('load ',r_dir,'owfields_',num2str(ti+(j-1)*Plang),'.mat');
        load_dle=strcat('load ',r_dir,'dlefields_',num2str(ti+(j-1)*Plang),'.mat');
        load_meso=strcat('load ',r_dir,'mesofields_',num2str(ti+(j-1)*Plang),'.mat');
        load_geo=strcat('load ',r_dir,'geofields_',num2str(ti+(j-1)*Plang),'.mat');
        load_fte=strcat('load ',r_dir,'ftefields_',num2str(ti+(j-1)*Plang),'.mat');
        load_ai=strcat('load ',r_dir,'aifields_',num2str(ti+(j-1)*Plang),'.mat');
        load_ftcs=strcat('load ',r_dir,'ftcsfields_',num2str(ti+(j-1)*Plang),'.mat');
        
        eval(load_data)
        xnow(1:Mlen(1),:)=xb{1};           
        ynow(1:Mlen(1),:)=yb{1};
        a=Mlen(1);
        if(size(Mlen,2)~=1)          
            for ii=2:size(Mlen,2)
                xnow(a+1:a+Mlen(ii),:)=xb{ii};
                ynow(a+1:a+Mlen(ii),:)=yb{ii};
                a=a+Mlen(ii);
            end
        end 
        for ti=1:Plang
            XX=xnow(:,ti);
            YY=ynow(:,ti);
            
            x0=1+round((((Nx-1)/diff(domain(1,:)))*(XX-domain(1,1))));
            y0=1+round((((Ny-1)/diff(domain(2,:)))*(YY-domain(2,1))));
            B0=Nx*(y0-1)+x0;
    
            if(strcmp(res_type,'FTLE'))
                eval(load_dle)
    
                [trac_val]=interpn(xs,ys,dle,XX,YY); % Interpolate to the grid values and dtermine region
                ind1=find(trac_val<=(mean(dle(:))-std(dle(:))));
                ind2=find(trac_val>=(mean(dle(:))+std(dle(:))));
                ind3=find(abs(trac_val-mean(dle(:)))<std(dle(:)));
                trac_val(ind1)=-1;
                trac_val(ind2)=1;
                trac_val(ind3)=0;
        
            elseif(strcmp(res_type,'OW'))
                eval(load_ow)
    
                [trac_val]=interpn(xs,ys,Q,XX,YY); % Interpolate to the grid values and dtermine region
                ind1=find(trac_val<=(mean(Q(:))-std(Q(:))));
                ind2=find(trac_val>=(mean(Q(:))+std(Q(:))));
                ind3=find(abs(trac_val-mean(Q(:)))<std(Q(:)));
                trac_val(ind1)=-1;
                trac_val(ind2)=1;
                trac_val(ind3)=0;
        
            elseif(strcmp(res_type,'MESO'))
                eval(load_meso)
                
                [trac_val]=interpn(xs,ys,meso,XX,YY); % Interpolate to the grid values and dtermine region
                % Put the correct threshold for MESO here.
                % I really need to edit this...
                ind1=find(trac_val<=(mean(dle(:))-std(dle(:))));
                ind2=find(trac_val>=(mean(dle(:))+std(dle(:))));
                ind3=find(abs(trac_val-mean(dle(:)))<std(dle(:)));
                trac_val(ind1)=-1;
                trac_val(ind2)=1;
                trac_val(ind3)=0;
        
            elseif(strcmp(res_type,'FTE'))
                eval(load_fte)
    
                [trac_val]=interpn(xs,ys,FTE,XX,YY); % Interpolate to the grid values and dtermine region
                ind1=find(trac_val<=(mean(FTE(:))-std(FTE(:))));
                ind2=find(trac_val>=(mean(FTE(:))+std(FTE(:))));
                ind3=find(abs(trac_val-mean(FTE(:)))<std(FTE(:)));
                trac_val(ind1)=-1;
                trac_val(ind2)=1;
                trac_val(ind3)=0;
        
            elseif(strcmp(res_type,'GEO'))     
                eval(load_geo);
                trac_val=GEO(B0(:));
        
            elseif(strcmp(res_type,'FTCS'))     
                eval(load_ftcs)
                trac_val=FTCS(B0(:));

            elseif(strcmp(res_type,'AI'))
                eval(load_ai)
                trac_val=AI(B0(:));
            else
                eval(load_ps);
                trac_val=PS(B0(:));
            end
            
            xl=[xl,trac_val];
            ti
        end
        
        if(j~=Narteng)
            for i=1:Ny
                xi=xl(1+(i-1)*(K*K*Nx):(K*K*Nx)*i,:);
                saver=strcat('save ../tmp/temp',num2str(j),'_',num2str(i),'.mat xi');
                eval(saver);
            end
            clear xl yl
        else
            for i=Ny:-1:1
                xz=[];
                i
                for ii=1:Narteng-1
                    loader=strcat('load ../tmp/temp',num2str(ii),'_',num2str(i),'.mat');
                    eval(loader);
                    xz=[xz,xi];
                end
                xz=[xz,xl(1+(i-1)*(K*K*Nx):(K*K*Nx)*i,:)];
                saver=strcat('save ',s_dir,noise_type,'/',res_type,'_zone',num2str(i),'.mat xz')
                eval(saver);
            end
        end
    end

    case 'moms'
        save_full=strcat('save ',s_dir,noise_type,'/full_moments.mat', ' meanx meany Sxx Syy Sxy Syx');
        
        tic
        for i=1:Nx*Ny
            xa=[]; ya=[];
            savefile=strcat('save ../tmp/temp_cell_',num2str(i),'.mat', ' xa ya');
            eval(savefile);
        end

        for i=1:Narteng
            xnow=[]; ynow=[];
            load_data=strcat('load ',r_dir,'final_positions',noise_type,num2str(i),'.mat');
            eval(load_data);
    
            xnow(1:Mlen(1),:)=xb{1};              
            ynow(1:Mlen(1),:)=yb{1};
            a=Mlen(1);
            if(size(Mlen,2)~=1)           
                for ii=2:size(Mlen,2)
                    xnow(a+1:a+Mlen(ii),:)=xb{ii};
                    ynow(a+1:a+Mlen(ii),:)=yb{ii};
                    a=a+Mlen(ii);
                end
            end 
            if(i==1)
                Xr=1+round(((Nx-1)/diff(domain(1,:)))*(xnow(:,1)-domain(1,1)));
                Yr=1+round(((Ny-1)/diff(domain(2,:)))*(ynow(:,1)-domain(2,1))); 
                Br=(Nx)*(Yr-1)+(Xr);
                BR=unique(Br);
                for j=1:length(BR)
                   index{j}=find(Br==BR(j));
                end 
            end
    
            for j=1:length(BR)
                loader=strcat('load ../tmp/temp_cell_',num2str(BR(j)),'.mat');
                eval(loader);
        
                xa=[xa,xnow(index{j},:)];
                ya=[ya,ynow(index{j},:)];
                saver=strcat('save ../tmp/temp_cell_',num2str(BR(j)),'.mat', ' xa ya');
                eval(saver);                
            end
            %clear xa ya
            i
        end
disp('done...and...');
toc
        tic
        meanx=zeros(Nx*Ny,n_cut*Narteng);
        meany=meanx;
        Sxx=meanx;
        Syy=meanx;
        Sxy=meanx;
        Syx=meanx;    
        for i=1:Nx*Ny
            loader=strcat('load ../tmp/temp_cell_',num2str(i),'.mat');
            eval(loader);
            i
            for j=1:n_cut*Narteng
                meanx(i,j)=mean(xa(:,j));
                meany(i,j)=mean(ya(:,j));
                C=cov(xa(:,j),ya(:,j));
                Sxx(i,j)=C(1,1);
                Sxy(i,j)=C(1,2);
                Syx(i,j)=C(2,1);
                Syy(i,j)=C(2,2);    
            end
        end
        eval(save_full)
        
        for j=1:n_cut*Narteng
            j
            mx=meanx(:,j);
            my=meany(:,j);
            sxx=Sxx(:,j);
            sxy=Sxy(:,j);
            syx=Syx(:,j);
            syy=Syy(:,j);
 
            save_part=strcat('save ',s_dir,noise_type,'/moms_',num2str(j),'.mat',' mx my sxx sxy syx syy');
            eval(save_part)
        end
        toc
end
toc