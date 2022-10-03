%%

load FullCode/prob_sets.mat

load FullCode/in_progress/BickFullAI_sets.mat

for t0=1:20
ind1=[]; %Would be able to then do a loop, and sort this by number of sets and sort(unique)
ind1{1}=find(AI(:,:,t0)==2);
ind1{2}=find(AI(:,:,t0)==-1);
ind1{3}=find(AI(:,:,t0)==1);
ind1{4}=find(AI(:,:,t0)==0);

ind2=[];
ind2{1}=find(FTCS(:,:,t0)==0);
ind2{2}=find(FTCS(:,:,t0)==1);
ind2{3}=find(FTCS(:,:,t0)==2);

Acheck=ones(144,64);
%Elliptic
for i=1:length(ind1{1})
    ind=find(ind1{1}(i)==ind2{2});
    if(~isempty(ind))
        Acheck(ind2{2}(ind))=0;
    end
end

for i=1:length(ind1{2})
    ind=find(ind1{2}(i)==ind2{3});
    if(~isempty(ind))
        Acheck(ind2{3}(ind))=0;
    end
end

for i=1:length(ind1{3})
    ind=find(ind1{3}(i)==ind2{1});
    if(~isempty(ind))
        Acheck(ind2{1}(ind))=0;
    end
end

for i=1:length(ind1{4})
    ind=find(ind1{4}(i)==ind2{3});
    if(~isempty(ind))
        Acheck(ind2{3}(ind))=0;
    end
end

subplot(3,1,1)
pcolor(FTCS(:,:,t0)');shading flat;
subplot(3,1,2)
pcolor(AI(:,:,t0)');shading flat;
subplot(3,1,3)
pcolor(Acheck');shading flat;
drawnow 
shg






end
%%


save Bickley_prob_parts.mat AI FTCS
