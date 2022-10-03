function eq_adv(XS,YS,varargin)
global domain resolution Nx Ny dxx
global s_dir
global comp_type noise_type pu_type
global Period INTTIME CaseSTART N_Basetime dt pt nsteps n_cut
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


%for tt=1:length(e_vec)
for ti=1:N_Basetime
    start_time=(ti-1)*pt+CaseSTART;
    ti
    %tt=1;
    %e1=e_vec(tt); %Can be removed if e1 is fixed
    xnow = XS;
    ynow = YS;
            
    for timeunit=1:INTTIME/dt;
        x0=xnow;y0=ynow;t0=start_time+(timeunit-1)*dt;
        if(strcmp(pu_type,'GPU'))
            [u1,v1]=arrayfun(@derivative,t0,x0,y0);
            [t1,x1,y1]=arrayfun(@steps,x0,y0,t0,u1,v1,.5,dt);
                    
            [u2,v2]=arrayfun(@derivative,t1,x1,y1);
            [t2,x2,y2]=arrayfun(@steps,x0,y0,t0,u2,v2,.5,dt);
                   
            [u3,v3]=arrayfun(@derivative,t2,x2,y2);
            [t3,x3,y3]=arrayfun(@steps,x0,y0,t0,u3,v3,1,dt);

            [u4,v4]=arrayfun(@derivative,t3,x3,y3);
            [xnow,ynow]=arrayfun(@last_step,xnow,ynow,u1,u2,u3,u4,v1,v2,v3,v4,dt);  
        else
            [u1,v1]=derivative(t0,x0,y0);
            [t1,x1,y1]=steps(x0,y0,t0,u1,v1,.5,dt);
                    
            [u2,v2]=derivative(t1,x1,y1);
            [t2,x2,y2]=steps(x0,y0,t0,u2,v2,.5,dt);
                   
            [u3,v3]=derivative(t2,x2,y2);
            [t3,x3,y3]=steps(x0,y0,t0,u3,v3,1,dt);

            [u4,v4]=derivative(t3,x3,y3);
            [xnow,ynow]=last_step(xnow,ynow,u1,u2,u3,u4,v1,v2,v3,v4,dt);
        end
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
                flag1=ceil(flag/nsteps);
                flag1=mod(flag1,n_cut);
                xz(:,flag1+1)=xtemp(:);        % Tracer array at field times
                yz(:,flag1+1)=ytemp(:);
                if(mod(flag1+1,n_cut)==0)
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
        dle_2d;
        meso_2d;
    end
end           
end

function [tt,xx,yy]=steps(x,y,t,u,v,b,dt)
    xx=x+u*dt*b;
    yy=y+v*dt*b;
    tt=t+b*dt;
end

function [xnow,ynow]=last_step(xnow,ynow,u1,u2,u3,u4,v1,v2,v3,v4,dt)
    xnow=xnow+dt*(u1+2*u2+2*u3+u4)/6;
    ynow=ynow+dt*(v1+2*v2+2*v3+v4)/6;
end
