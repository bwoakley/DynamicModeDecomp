close all
clear all
clc
addpath ../mytoolbox/
r_dir=strcat(' data/'); %read directory
s_dir=strcat(' data/'); %save directory
cut=1;

loadfile=strcat(r_dir,'initial_positions1.mat');
sss=strcat('load',loadfile);
eval(sss);
loadfile=strcat(r_dir,'final_positions1.mat');
sss=strcat('load',loadfile);
eval(sss);

    num_cuts=15;

    
for i=1:size(xa,2)
    x0=xa{i};
    y0=ya{i};
    xf=xb{i};
    yf=yb{i};
    
    if(i==1)
        s=floor(size(x0,1)/num_cuts);
        for j=1:num_cuts-1;
            arc0(j)=arclength(x0(1+s*(j-1):s*j),y0(1+s*(j-1):s*j));
            arcf(j)=arclength(xf(1+s*(j-1):s*j),yf(1+s*(j-1):s*j));
        end
        arc0(num_cuts)=arclength(x0(s*(num_cuts-1)+1:end),y0(s*(num_cuts-1)+1:end));
        arcf(num_cuts)=arclength(xf(s*(num_cuts-1)+1:end),yf(s*(num_cuts-1)+1:end));
        
    else
        figure(i)
        plot(xa{1},ya{1},'k'); hold on
        % Search for the most similar point on the geodesic curves, for the
        % endpoints of the segments used in the Poincare section.
        for j=1:num_cuts-1
            idx1=knnsearch([x0,y0],[xa{1}(1+s*(j-1)),ya{1}(1+s*(j-1))]);
            idx2=knnsearch([x0,y0],[xa{1}(s*j),ya{1}(s*j)]);
            
            if(j==cut)
            plot(x0(idx1),y0(idx1),'gp')
            plot(xf(idx1),yf(idx1),'rp')
            end
            if(idx1<idx2)
                arc0(j)=arclength(x0(idx1:idx2),y0(idx1:idx2));
                arcf(j)=arclength(xf(idx1:idx2),yf(idx1:idx2));
                if(j==cut)
                %plot(x0(idx1:idx2),y0(idx1:idx2),'g')
                plot(x0(:),y0(:),'g')
                plot(xf(:),yf(:),'r')
                %plot(xf(idx1:idx2),yf(idx1:idx2),'r')
                end
            else
                arc0(j)=arclength(x0([idx1:end,1:idx2]),y0([idx1:end,1:idx2]));
                arcf(j)=arclength(xf([idx1:end,1:idx2]),yf([idx1:end,1:idx2]));
                if(j==cut)
                %plot(x0([idx1:end,1:idx2]),y0([idx1:end,1:idx2]),'g')
                plot(x0(:),y0(:),'g')
                plot(xf(:),yf(:),'r')
                %plot(xf([idx1:end,1:idx2]),yf([idx1:end,1:idx2]),'r')
                end
            end
        end
        idx1=knnsearch([x0,y0],[xa{1}(s*(num_cuts-1)+1),ya{1}(s*(num_cuts-1)+1)]);
        idx2=knnsearch([x0,y0],[xa{1}(end),ya{1}(end)]);
        if(idx1<idx2)
            arc0(num_cuts)=arclength(x0(idx1:idx2),y0(idx1:idx2));
            arcf(num_cuts)=arclength(xf(idx1:idx2),yf(idx1:idx2));
            if(j==cut)
            %plot(x0(idx1:idx2),y0(idx1:idx2),'g')
                plot(x0(:),y0(:),'g')
                plot(xf(:),yf(:),'r')
                %plot(xf(idx1:idx2),yf(idx1:idx2),'r')
            end
        else
            arc0(num_cuts)=arclength(x0([idx1:end,1:idx2]),y0([idx1:end,1:idx2]));
            arcf(num_cuts)=arclength(xf([idx1:end,1:idx2]),yf([idx1:end,1:idx2]));
            if(j==cut)
            %plot(x0([idx1:end,1:idx2]),y0([idx1:end,1:idx2]),'g')
                plot(x0(:),y0(:),'g')
                plot(xf(:),yf(:),'r')
                %plot(xf([idx1:end,1:idx2]),yf([idx1:end,1:idx2]),'r')
            end
        end
    end
    arc0(num_cuts+1)=arclength(x0,y0);
    arcf(num_cuts+1)=arclength(xf,yf);    
%Now the arc lengths of the smaller line segments have been determined,
%and the difference between before and after being advected are listed.
    lambda_guess1{i}=(arcf-arc0)./arcf;
    lambda_guess2{i}=(arcf-arc0)./arc0;
    
%     figure(i)
%     plot(x0,y0,'k'); hold on
%     plot(xf,yf,'r')
end

savefile=strcat(s_dir,'ratios1.mat');
sss=strcat('save',savefile,' lambda_guess1 lambda_guess2');
eval(sss);

break
%% Compare to the other things.

N=30; %For e1=.002;
%N=60; %For e1=upper;
lamb_vec=linspace(.85,1.15,N);

loadfile=strcat('load',r_dir,'ratios1.mat');
eval(loadfile);


    
c vvvvvvvvvvvvvvvvvvvvvmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmvvvvvvvvvvvvvvvvvv