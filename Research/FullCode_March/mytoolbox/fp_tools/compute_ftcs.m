function [B11,B12,C11,C12,youtmax]=compute_ftcs(Di2,A2,C1,B1,nsteps)

q=(Di2'*A2)';
ind=find(q==0);
if(~isempty(ind))
    index=true(1,size(A2,2));
    index(ind)=false;
    A2=A2(:,index);
end
q=(Di2'*A2)';
ScaleLeft=diag((Di2.^(1/2)));
ScaleRight=diag((q.^(1/2)));
L=ScaleLeft*A2/(ScaleRight);
[U,~,V]=svds(L,2);
leftvecs=(U'/(ScaleLeft))';
rightvecs=(V'/(ScaleRight))';
[youtmax,~,~,~,~,lvec0,rvec0]=threshold_coherent_nonsquareP(A2,leftvecs(:,2),rightvecs(:,2),nsteps,Di2);

B11=B1((lvec0==1));
B12=B1((lvec0==2));
C11=C1((rvec0==1));
C12=C1((rvec0==2));

end


function [youtmax,ioutmax,iout1,iout2,iout2match,leftthreshvec,rightthreshvec]=threshold_coherent_nonsquareP(P,leftvec,rightvec,step,p)
%P is the transition matrix, p is the invariant measure, normalised to sum to 1.
%leftvec, rightvec come from coherent_vectors.m
%step is the integer number of boxes that are stepped through in the line optimisation (an integer;  1 is accurate (steps through every single box) and slower, larger
%numbers are faster and less accurate). 
%assumes P is row stochastic
if size(p,1)>size(p,2),
    p=p';
end
p=p/sum(p);
q=p*P;
n=size(P,1);
n2=size(P,2);

%% descend
[~,id1]=sort(leftvec,'descend');
[~,id2]=sort(rightvec,'descend');
pdesc=p(id1);
Pdesc=P(id1,id2);
ratiodesc=zeros(1,n);
pdescsum=cumsum(p(id1));
qdescsum=cumsum(q(id2));
[~,idmax]=min(abs(pdescsum-0.5));

id2match=zeros(size(id1));
for i=1:step:idmax,
    [~,jmin]=min(abs(pdescsum(i)-qdescsum));
    id2match(i)=jmin; 
    mass_retain=sum(pdesc(1:i)*Pdesc(1:i,1:jmin));
    mass_total=pdescsum(i);
    ratiodesc(i)=mass_retain/mass_total;
end

%% ascend
[~,ia1]=sort(leftvec,'ascend');
[~,ia2]=sort(rightvec,'ascend');
pasc=p(ia1);
Pasc=P(ia1,ia2);
ratioasc=zeros(1,n);
pascsum=cumsum(p(ia1));
qascsum=cumsum(q(ia2));
[~,iamax]=min(abs(pascsum-0.5));

ia2match=zeros(size(id1));
for i=1:step:iamax,
    [~,jmin]=min(abs(pascsum(i)-qascsum));
    ia2match(i)=jmin; 
    mass_retain=sum(pasc(1:i)*Pasc(1:i,1:jmin));
    mass_total=pascsum(i);
    ratioasc(i)=mass_retain/mass_total;
end

%% plotting
[ydesc,idesc]=max(ratiodesc);
[yasc,iasc]=max(ratioasc);

%figure
%hold on;
if ydesc>yasc,
    youtmax=ydesc;
    ioutmax=idesc;
    iout1=id1;
    iout2=id2;
    iout2match=id2match;
%    plot(ioutmax,youtmax,'k*','markersize',10)
else
    youtmax=yasc;
    ioutmax=iasc;
    iout1=ia1;
    iout2=ia2;
    iout2match=ia2match;
%    plot(n-ioutmax+1,youtmax,'k*','markersize',10)
end
%plot(1:n,ratiodesc,'b.');hold on; plot(n:-1:1,ratioasc,'r.')

leftthreshvec=zeros(n,1);
rightthreshvec=zeros(n2,1);
leftthreshvec(iout1(1:ioutmax))=1;
leftthreshvec(setdiff(1:n,iout1(1:ioutmax)))=2;
rightthreshvec(iout2(1:iout2match(ioutmax)))=1;
rightthreshvec(setdiff(1:n2,iout2(1:iout2match(ioutmax))))=2;
end