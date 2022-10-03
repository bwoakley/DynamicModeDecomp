
r_dir=strcat(' ../data_revs/'); %read tracer directory
s_dir=strcat(' ../data_revs/'); %save tracer directory

runtype='AI'; %Or 'FTCS'
ti=3;
tf=20;


if(tf>ti)
    dir=1;
else
    dir=-1;
end
for t0=ti:dir:tf
    tic
    loader1 = strcat('load ',r_dir,'FP_set_',num2str(t0),'.mat Ab Di1 Nx Ny');
    loader2 = strcat('load ',r_dir,'FP_set_',num2str(t0),'.mat B0 Br');
    
    save_ai = strcat('save -v7.3',s_dir,'FP_ai_',num2str(t0),'.mat',' I Ny Nx');
    save_ftcs = strcat('save -v7.3 ',s_dir,'FP_ftcs_',num2str(t0),'.mat A');

    eval(loader1)
    Ab=sparse(Ab);
    switch runtype
        case 'AI'            
            
        [Uxx,Vxx]=eigs(Ab',3);
        p=Uxx(:,1); PI = p;
        p=Di1./sum(Di1);
        Phat = spdiags((1./PI),0,length(PI),length(PI))*Ab'*spdiags(PI,0,length(PI),length(PI));
        R=((Ab+Phat)/2);
        clear Ab Phat

        q=5;  %I should define these outside of the loop
        s=Nx*Ny;
        m=2;
        [ur,vr]=eigs(R,q+1);
        clear R
        [vir,id]=sort(diag(vr),'descend'); ur=ur(:,id);

        c=17;
        q=2;
        x=ur(:,2:q+1);
        for i=1:size(x,2)
            x(:,i)=(x(:,i))./norm(x(:,i));
        end
        % call picard;
        [U,v_his]=picard_its(x',c,s,m,p);
        [~,I] = max(U); 
        flag=1;
        while(flag<5)
            [~,I] = max(U);
            for ii=1:c
                ind=find(I==ii);
                rho_=numel(ind);
                if(rho_/(s)<(1/(c+5)))
                    U(ii,:)=0;
                end
            end
            flag=flag+1;
        end
        [~,I] = max(U);
        eval(save_ai)
        
    case 'FTCS'
        lvl=0; rho0=.955; 
    
        B1_parts{1}=B;
        C1_parts{1}=B;
        keeper=1;
        Keep=[];
        while(~isempty(B1_parts) && lvl<5)
            flag=0; B2_parts=[]; C2_parts=[];
            for jj=1:numel(B1_parts)
                R=B1_parts{jj};
                R2=C1_parts{jj};
        
                if(numel(R)<10 || numel(R2)<10)
                    %skip and store the parent partition.
                    rho=.4;
                else
                    if(lvl==0)
                        %skip the reforming of Ab and Di1
                        A2=Ab;
                        Di2=Di1;
                    else
                        index=false(1,size(Ab,1));
                        index(R)=true;
                        Di2=Di1(index,:);
                        A2=Ab(index,:);
                        index=false(1,size(Ab,2));
                        index(R2)=true;
                        A2=A2(:,index);    
                    end
                
                    if(lvl<2)
                        A2=sparse(A2);
                        nsteps=1;
                    else
                        nsteps=1;
                    end
                    [B11,B12,C11,C12,rho]=compute_ftcs(Di2,A2,R2,R,nsteps);
                end            
                rho
                if(rho>=rho0)
                    B2_parts{flag+1}=B11;
                    B2_parts{flag+2}=B12;
                    C2_parts{flag+1}=C11;
                    C2_parts{flag+2}=C12;
                    flag=flag+2;
                else
                    Keep{keeper}=R;
                    keeper=keeper+1;
                end
            end
            B1_parts=[];
            C1_parts=[];
            B1_parts=B2_parts;
            C1_parts=C2_parts;
            lvl=lvl+1;
        end

        A=zeros(Ny,Nx);
        for i=1:numel(B2_parts)
            A(B2_parts{i})=i;
        end
        for j=1:numel(Keep)
            A(Keep{j})=numel(B2_parts)+j;
        end
        figure;pcolor(A);shading flat;shg
    
        eval(save_ftcs)
    end
    toc
end
