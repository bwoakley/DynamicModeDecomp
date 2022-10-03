% This is a filament check for flow field partitioning metrics.

%Need to figure out how to apply masks
global dxx Nx Ny domain resolution
global s_dir rd_dir
global inter
global comp_type noise_type pu_type
global Period INTTIME CaseSTART N_Basetime dt pt nsteps n_cut ft
global periodic

r_dir='';
load globs.mat

prompt='(GPU) or (CPU)? ';
pu_type=input(prompt);

prompt='What start time? ';
t0=input(prompt);

prompt='Integrate for how many periods? ';
nP=input(prompt);

INTTIME=20*nP;
periodic=[1 1];


file_ind=ceil(t0/n_cut);
flag=mod(t0,n_cut);
if(flag==0)
   file_num=file_ind;
   flag=n_cut;
else
   file_num+file_in+1;
end

%This should be changeable...
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

xnow=xz(:,flag);
ynow=yz(:,flag);
if(periodic(1))
    xnow=mod(xnow-domain(1,1),diff(domain(1,:)))+domain(1,1);
end
if(periodic(2))
    ynow=mod(ynow-domain(2,1),diff(domain(2,:)))+domain(2,1);
else
    ynow(ynow>domain(2,2))=domain(2,2);
    ynow(ynow<domain(2,1))=domain(2,1);
end

%xnow and ynow are the original locations... to be used as input for the advection.
flag=flag+.3; %This will ensure that no data will be overwritten because the checks wont succeed.
N_Basetime=1;
CaseSTART=(t0-1)*(pt/Period);
nsteps=100000;
%At this point, will have 
% xnow ynow xx yy flag file_num periodic INTTIME dt(globs).
switch eq_type
    case 'data'
    	data_adv(xnow,ynow,xx,yy); %Use data to compute traj.
    
    case 'eq'
    	eq_adv(xnow,ynow); %Use equations to compute traj.
end
%Then save the output...



%This bit of post-processing is gomnma m,e





%Possibly do some plotting


