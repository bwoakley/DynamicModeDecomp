%% This is purely to add buffer regions to already partitioned fields.

%Specify the region value (reg_val), and the buffer value (buf_val).
reg_val=1;
buf_val=40;

%Specify the read and write (r_dir and s_dir, resp.)
r_dir=strcat(' ../somewhere/else/');
s_dir=strcat(' ../../final_fields/');

field_type='ftle'
ti=1;
tf=400;
dir=1;

% Need to specify Nx and Ny somehow...(read-in globs.mat)

for t0=ti:dir:tf
    loader=strcat('load ',r_dir,field_type,'_',num2str(t0),'.mat');
    saver=strcat('save ',s_dir,field_type,'_',num2str(t0),'.mat field');
    eval(loader);
    A=field;
    
    B=A;
    for i=1:Nx
        for j=1:Ny
            if(A(i,j)==reg_val)
                if(j<Ny) 
                    if(A(i,j+1)~=reg_val)
                        B(i,j+1)=buf_val;
                    end
                end
                if(j>1)
                    if(A(i,j-1)~=reg_val) 
                        B(i,j-1)=buf_val;
                    end
                end
                if(i<Nx)
                    if(A(i+1,j)~=reg_val) 
                        B(i+1,j)=buf_val;
                    end
                end
                if(i>1)
                    if(A(i-1,j)~=reg_val) 
                        B(i-1,j)=buf_val;
                    end
                end
            end
        end
    end
    field=B;
    eval(saver);
end