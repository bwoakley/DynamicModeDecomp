function [youtmax,ioutmax,iout1,iout2,iout2match,leftthreshvec]=threshold_invariant(leftvec,step,p)

%P is the transition matrix, p is the invariant measure, normalised to sum
%to 1.
%leftvec, rightvec come from coherent_vectors.m
%step is the integer number of boxes that are stepped through in the line optimisation (an integer;  1 is accurate (steps through every single box) and slower, larger
%numbers are faster and less accurate). 
%assumes P is row stochastic

if size(p,1)>size(p,2),
    p=p';
end
p=p/sum(p);

%% descend

[cout1,id1]=sort(leftvec,'descend');
pdesc=p(id1);
ratiodesc=zeros(1,n);
pdescsum=cumsum(p(id1));

[dum,idmax]=min(abs(pdescsum-0.5));

id2match=zeros(size(id1));
for i=1:step:idmax,
    [dum,jmin]=min(abs(pdescsum(i)-qdescsum));
    id2match(i)=jmin; 
    mass_retain=sum(pdesc(1:i)*Pdesc(1:i,1:jmin));
    mass_total=pdescsum(i);
    ratiodesc(i)=mass_retain/mass_total;
end

%% ascend

[cout1,ia1]=sort(leftvec,'ascend');
pasc=p(ia1);
ratioasc=zeros(1,n);
pascsum=cumsum(p(ia1));

[dum,iamax]=min(abs(pascsum-0.5));

ia2match=zeros(size(id1));
for i=1:step:iamax,
    [dum,jmin]=min(abs(pascsum(i)-qascsum));
    ia2match(i)=jmin; 
    mass_retain=sum(pasc(1:i)*Pasc(1:i,1:jmin));
    mass_total=pascsum(i);
    ratioasc(i)=mass_retain/mass_total;
end

%% plotting

[ydesc,idesc]=max(ratiodesc);
[yasc,iasc]=max(ratioasc);

figure
hold on;

if ydesc>yasc,
    youtmax=ydesc;
    ioutmax=idesc;
    iout1=id1;
    iout2=id2;
    iout2match=id2match;
    plot(ioutmax,youtmax,'k*','markersize',10)
else
    youtmax=yasc;
    ioutmax=iasc;
    iout1=ia1;
    iout2=ia2;
    iout2match=ia2match;
    plot(n-ioutmax+1,youtmax,'k*','markersize',10)
end

plot(1:n,ratiodesc,'b.');hold on; plot(n:-1:1,ratioasc,'r.')


    
leftthreshvec=zeros(n,1);
rightthreshvec=zeros(n2,1);
leftthreshvec(iout1(1:ioutmax))=1;
leftthreshvec(setdiff(1:n,iout1(1:ioutmax)))=2;
rightthreshvec(iout2(1:iout2match(ioutmax)))=1;
rightthreshvec(setdiff(1:n2,iout2(1:iout2match(ioutmax))))=2;