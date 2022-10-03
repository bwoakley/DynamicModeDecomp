%This code uses the MinBoundSuite
clear all; close all
addpath ../mytoolbox/MinBoundSuite/

r_dir=strcat(' data_refine/'); %read directory
s_dir=strcat(' data_refine/'); %save directory

forwardLcsColor = 'r';
backwardLcsColor = 'b';
shearLcsColor = [0,.6,0];
strainLcsColor = [0,0,0];
sLcsColor = [0,0,1];

%First_time=1; lets you see the geodesics for the lambda layer
%First_time=2; picks the best layers and composes them.
first_time=1;

prompt1='What ti?';
ti=input(prompt1);
domain=[[-1,1]*2.5*pi;[-1,1]*3];
tar1=1; 
tar2=1;
for lam_i=61:120
    loadfile = strcat(r_dir,'geo_',num2str(ti),'_lambda',num2str(lam_i),'.mat');
    ssss = strcat('load',loadfile);
    eval(ssss)
    
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
end

if(first_time==1)
    lam_vec=[1:60]+60;
    s_vec(1)=0; s_vec(2)=0;
    narf=0;
else
    prompt= 'what lambda values?';
    lam_vec=input(prompt);
    prompt= 'Should we skip some?';
    s_vec=input(prompt);
    prompt= 'What lambda layer is it in?';
    narf=input(prompt);
end


for lam_j=1:size(lam_vec,2)
    loadfile = strcat(r_dir,'geo_',num2str(ti),'_lambda',num2str(lam_vec(lam_j)),'.mat');
    ssss = strcat('load',loadfile);
    eval(ssss)

    figure(lam_j);clf;
    pcolor(linspace(domain(1,1),domain(1,2),144*4*2),linspace(domain(2,1),domain(2,2),64*4*2),dle); hold on; shading interp
    xlim([domain(1,1)-diff(domain(1,:)) domain(1,2)+diff(domain(1,:))])
    colormap(flipud(gray))
    drawnow
    
    for j = 1:nPoincareSection
        hClosedOrbitsEtaPos = arrayfun(@(i)plot(closedOrbits{i}{1}{end}(:,1),closedOrbits{i}{1}{end}(:,2)),1:size(closedOrbits,2));
        set(hClosedOrbitsEtaPos,'color',shearLcsColor)
        set(hClosedOrbitsEtaPos,'linewidth',2)
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
        
        hClosedOrbitsEtaNeg = arrayfun(@(i)plot(closedOrbits{i}{2}{end}(:,1),closedOrbits{i}{2}{end}(:,2)),1:size(closedOrbits,2));
        set(hClosedOrbitsEtaNeg,'color',strainLcsColor)
        set(hClosedOrbitsEtaNeg,'linewidth',2)
        
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
end

if(first_time~=1)   
    real_keep1=zeros(nPoincareSection,1);
    for j = 1:nPoincareSection
        fixed_list=[];
        same_list=[];
        center_list=[];
        n=tar1(j);
        if(n==1)%There is no geodesic for this i.c.
        elseif(n==2)%This is the maximum geodesic for that i.c.
            real_keep1(j,:)=1;
        else
            for ii=1:n-2
                X1=n_closed{j}{ii}(:,1);
                Y1=n_closed{j}{ii}(:,2);
                [c1,r1]=minboundcircle(X1,Y1);
                for jj=ii+1:n-1
                    X2=n_closed{j}{jj}(:,1);
                    Y2=n_closed{j}{jj}(:,2);
                    [c2,r2]=minboundcircle(X2,Y2);
                    radi=max(r2,r1);
                    if(r2>=r1)
                        fixed_list=[fixed_list;jj,ii,radi];
                    else
                        fixed_list=[fixed_list;ii,jj,radi];
                    end
                end
            end
            if(~isempty(fixed_list))
                check=unique(fixed_list(:,2));
                keep=[];
                for i=1:length(check)
                    ind=find(check(i)==fixed_list(:,2));
                    [~,v]=max(fixed_list(ind,3));
                    keep(i)=fixed_list(ind(v),1);
                end
                real_keep1(j,:)=unique(keep)';
            else
                real_keep1=[];
            end
        end
    end
        
    real_keep2=zeros(nPoincareSection,1);
    for j=1:nPoincareSection
        fixed_list=[];
        same_list=[];
        center_list=[];
        n=tar2(j);
        if(n==1)%There is no geodesic for this i.c.
        elseif(n==2)%This is the maximum geodesic for that i.c.
            real_keep2(j,:)=1;
        else
            for ii=1:n-2
                X1=p_closed{j}{ii}(:,1);
                Y1=p_closed{j}{ii}(:,2);
                [c1,r1]=minboundcircle(X1,Y1);
                for jj=ii+1:n-1
                    X2=p_closed{j}{jj}(:,1);
                    Y2=p_closed{j}{jj}(:,2);
                    [c2,r2]=minboundcircle(X2,Y2);
                    radi=max(r2,r1);
                    if(r2>=r1)
                        fixed_list=[fixed_list;jj,ii,radi];
                    else
                        fixed_list=[fixed_list;ii,jj,radi];
                    end
                end
            end
            if(~isempty(fixed_list))
                check=unique(fixed_list(:,2));
                keep=[];
                for i=1:length(check)
                    ind=find(check(i)==fixed_list(:,2));
                    [~,v]=max(fixed_list(ind,3));
                    keep(i)=fixed_list(ind(v),1);
                end
                real_keep2(j,:)=unique(keep)';
            else
                real_keep2=[];
            end
        end
    end
    for i=1:nPoincareSection
        geodesX{i}=[NaN,NaN];
        geodesY{i}=[NaN,NaN];
    end
    for i=1:nPoincareSection
        if(real_keep1(i)~=0)
            if(n_closed{i}{real_keep1(i)}(:,1)>=-9.25);
            geodesX{i}=n_closed{i}{real_keep1(i)};
            end
        end
        if(real_keep2(i)~=0)
            if(p_closed{i}{real_keep2(i)}(:,1)>=-9.25);
            geodesY{i}=p_closed{i}{real_keep2(i)};
            end
        end
    end
    
    hClosedOrbitsEtaPos = arrayfun(@(i)plot(geodesX{i}(:,1),geodesX{i}(:,2)),1:nPoincareSection);
    set(hClosedOrbitsEtaPos,'color',sLcsColor)
    set(hClosedOrbitsEtaPos,'linewidth',2)
    hClosedOrbitsEtaPos = arrayfun(@(i)plot(geodesY{i}(:,1),geodesY{i}(:,2)),1:nPoincareSection);
    set(hClosedOrbitsEtaPos,'color',sLcsColor)
    set(hClosedOrbitsEtaPos,'linewidth',2)

    MM=nPoincareSection;
    %savefile=strcat(s_dir,'geo_poss_chk',num2str(ti),'.mat');
    %ssss = strcat('save',savefile, ' geodesX geodesY MM domain resolution');
    %eval(ssss)
end