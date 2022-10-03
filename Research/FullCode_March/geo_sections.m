%% This code will read in the Poincare Section and the Geodesic information,
% and cut them into N segments. Save the cut up bits to be used for the
% advection code.

close all
clear all
clc
r_dir=strcat(' data/'); %read directory
s_dir=strcat(' data/'); %save directory

load ..\..\..\Downloads\PBND2\PoincBnd.mat
ti=1;
%xp=xpa; yp=ypa;
for lami=1:15
xxx=(mod(xp(:,ti:20:end)+2.5*pi,5*pi))-2.5*pi;
yyy=(yp(:,ti:20:end));

MM=6;
k = convhull(xxx(MM,:),yyy(MM,:));     
xps{lami}=xxx(MM,k);
yps{lami}=yyy(MM,k);
NN=floor(length(xps{lami})/10);

loadfile=strcat(r_dir,'geo_',num2str(ti),'_lambda',num2str(lami),'.mat');
sss=strcat(' load',loadfile);
eval(sss);
    
geodesX=[];
geodesY=[];
nPoincareSection=numel(closedOrbits);
list1=[];
list2=[];
tar1(1:nPoincareSection)=1;
tar2(1:nPoincareSection)=1;

for j = 1:nPoincareSection
    n_closed{j}=[];
    p_closed{j}=[];
end

lam_vec=[1:30];
    s_vec(1)=0; s_vec(2)=0;
    narf=0;
for j = 1:nPoincareSection
    if(~isnan(closedOrbits{j}{1}{end}))
        dd=diff(closedOrbits{j}{1}{end}(:,1));
        if(~isempty(find(dd>1e-6)))
            if(s_vec(1)==1 && s_vec(2)==j && narf==lam_vec(lam_j))
            else
                n_closed{j}{tar1(j)}=closedOrbits{j}{1}{end};
                tar1(j)=tar1(j)+1;
                list1=[list1,j];
            end
        end
    end 
    if(~isnan(closedOrbits{j}{2}{end}))
         dd=diff(closedOrbits{j}{2}{end}(:,1));
         if(~isempty(find(dd>1e-6)))
             if(s_vec(1)==2 && s_vec(2)==j && lam_vec(lam_j)==narf)
             else
                 p_closed{j}{tar2(j)}=closedOrbits{j}{2}{end};
                 tar2(j)=tar2(j)+1;
                 list2=[list2,j];
             end
         end
    end
end

MM=3;
k=convhull(p_closed{MM}{end}(:,1),p_closed{MM}{end}(:,2));
xps2{lami}=p_closed{MM}{end}(:,1);
yps2{lami}=p_closed{MM}{end}(:,2);


figure(lami)
plot(xps{lami},yps{lami},'.k','linewidth',2)
hold on
plot(xps2{lami},yps2{lami},'.r','linewidth',2)
end
 

break
%%
clear xa ya
keep=[1:11];
xa{1}=xps{1}(:); ya{1}=yps{1}(:);
for ii=1:size(keep,2)
    xa{1+ii}=xps2{keep(ii)};
    ya{1+ii}=yps2{keep(ii)};
end

savefile=strcat(s_dir,'initial_positions1.mat');
sss=strcat('save',savefile,' xa ya');
eval(sss);