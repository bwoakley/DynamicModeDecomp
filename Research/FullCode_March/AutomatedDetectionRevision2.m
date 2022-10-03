function [poincareSection,M]=AutomatedDetectionRevision2(poincareSection,cgStrain,cgEigenvector,showGraph,ftle)

xs=linspace(-2.5*pi,2.5*pi+5*pi,144*2+1); xs=xs(1:end-1);
ys=linspace(-3,3,65); ys=ys(1:end-1);
deltax=(xs(2)-xs(1));
[XS YS]=meshgrid(xs,ys);

showGraph='true'
if(showGraph)
figure(30)
imagesc(xs,ys,ftle)
colormap(flipud(gray))
hold on
end



%% STEP 1: Use an approximation in order to find the zero level curves...
deltai=[.65];%,.3];
for i=1:length(deltai)

ind=find(abs(cgStrain(1,1,:)-cgStrain(2,2,:)-cgStrain(1,2,:))<=deltai(i));
if(i==1)
    if(showGraph)
       plot(XS(ind),YS(ind),'gp')
    end
    larger_set=ind;
else
    if(showGraph)
       plot(XS(ind),YS(ind),'rp')
    end
    smaller_set=ind;
end


%% STEP 2: Use Delauny triangulation to set a minimum distance between
% singularities.

X=[]; Y=[];
for i=1:length(larger_set)
    X=[X, XS(larger_set(i))];
    Y=[Y, YS(larger_set(i))];
end
Z=[X',Y'];
[idx,dist]=knnsearch(Z,Z,'dist','euclidean','k',2);

keeper_set=[];
for i=1:length(larger_set)
    if(dist(i,2)>=1*deltax)
        keeper_set=[keeper_set,larger_set(i),larger_set(idx(i,2))];
    end
end

keeper_set=unique(keeper_set);


%% STEP 3: Determine the singularity type
r=.001; deltar=.001; N=1000;
Sing_type=zeros(length(keeper_set),N);
theta=linspace(0,2*pi,N);
for j=1:length(keeper_set)
    eta2_vec=cgEigenvector(keeper_set(j),3:4);
    r_vec=[r*cos(theta);r*sin(theta)];
    Sing_type(j,:)=abs(r_vec'*eta2_vec')./r;
    Sing_points_zero(j)=length(find(abs(Sing_type(j,:))<=deltar));
    Sing_points_one(j)=length(find(abs(Sing_type(j,:)-1)<=deltar));
    if(Sing_points_zero(j)==3 && Sing_points_one(j)==3)
        isTrisector(j)=1;
    else
        isTrisector(j)=0;
    end
end
end


%% STEP 4: Filter

%For now, just do the filtering on the larger set...use the KNN algorithm
X=[]; Y=[];
for i=1:length(keeper_set)
    X=[X, XS(keeper_set(i))];
    Y=[Y, YS(keeper_set(i))];
end
Z=[X',Y'];
[idx,dist]=knnsearch(Z,Z,'dist','euclidean','k',2);

% I need a more robust method for filtering the points...

% I should check the type of each point and make sure that they are
% compatible, and that if there are different with equal distances to the
% closest one, then I should look at the combination of those points....
keepers=[];

for i=1:length(keeper_set)
    if(isTrisector(i)+isTrisector(idx(i,2))<1 && dist(i,2)<=.5)
        keepers=[keepers,i,idx(i,2)];
    end
end


% Sort the keepers and remove doubles
keep=keepers;


X_keep=keeper_set(keep);
Y_keep=keeper_set(keep);

if(showGraph)
    plot(XS(X_keep),YS(Y_keep),'k+')
    hold on
    plot(Z(:,1),Z(:,2),'ro')
end


%% STEP 5: Use the mid-point of the remaining singularities to be the
% starting Poincare point
j=1;
X_loc=[]; Y_loc=[];
for i=1:2:length(X_keep)-1
    if(XS(X_keep(i))>9)
% %             X_loc(j)=(XS(X_keep(i+1))+XS(X_keep(i)))/2;
% %     Y_loc(j)=(YS(Y_keep(i+1))+YS(Y_keep(i)))/2;
% %     j=j+1;
    else
    X_loc(j)=(XS(X_keep(i+1))+XS(X_keep(i)))/2;
    Y_loc(j)=(YS(Y_keep(i+1))+YS(Y_keep(i)))/2;
    j=j+1;
;
    end
end

if(showGraph)
    plot(X_loc,Y_loc,'bp')
end


% This section defines the poincareSections, starting at the specified
% singularity locations, and end points pointing towards the zero.
for j=1:length(X_loc)
        poincareSection(j).endPosition = [X_loc(j), Y_loc(j); X_loc(j)-X_loc(j)/abs(X_loc(j)), Y_loc(j)];
        if(X_loc(j)>5 && X_loc(j)<6.3)
            poincareSection(j).endPosition = [X_loc(j), Y_loc(j); X_loc(j)+X_loc(j)/abs(X_loc(j)), Y_loc(j)];
        end
        
        if(X_loc(j)>7)
            if(X_loc(j)<8)
                poincareSection(j).endPosition = [X_loc(j), Y_loc(j); X_loc(j)-.5*X_loc(j)/abs(X_loc(j)), Y_loc(j)+.3];
            else
                poincareSection(j).endPosition = [X_loc(j), Y_loc(j); X_loc(j)+.5*X_loc(j)/abs(X_loc(j)), Y_loc(j)+.3];
            end
        end
end

M=length(X_loc);

