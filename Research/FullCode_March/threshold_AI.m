function [youtmax,ioutmax,iout1,threshvec]=threshold_AI(P,vec,step,p)

%P is the transition matrix, p is the invariant measure, normalised to sum
%to 1.
%vec is the vector you wish to threshold, eg. an eigenvector of P or R.
%step is the integer number of boxes that are stepped through in the line optimisation (an integer;  1 is accurate (steps through every single box) and slower, larger
%numbers are faster and less accurate). 
%assumes P is row stochastic

if size(p,1)>size(p,2),
    p=p';
end
p=p/sum(p);

N=length(vec);

%% descending

[cout1,id1]=sort(vec,'descend');

pdesc=p(id1);
Pdesc=P(id1,id1);
ratiodesc=zeros(1,N);
pdescsum=cumsum(pdesc);

[dum,idmax]=min(abs(pdescsum-0.5));

for i=1:step:idmax,
    mass_retain=sum(pdesc(1:i)*Pdesc(1:i,1:i));
    mass_total=pdescsum(i);
    ratiodesc(i)=mass_retain/mass_total;
end

%% ascending

[cout1,ia1]=sort(vec,'ascend');

pasc=p(ia1);
Pasc=P(ia1,ia1);
ratioasc=zeros(1,N);
pascsum=cumsum(pasc);

[dum,iamax]=min(abs(pascsum-0.5));

for i=1:step:iamax,
    mass_retain=sum(pasc(1:i)*Pasc(1:i,1:i));
    mass_total=pascsum(i);
    ratioasc(i)=mass_retain/mass_total;
end



%% plotting

[ydesc,idesc]=max(ratiodesc);
[yasc,iasc]=max(ratioasc);

figure;
plot(1:N,ratiodesc,'b.');hold on; plot(N:-1:1,ratioasc,'r.')

if ydesc>yasc,
    youtmax=ydesc;
    ioutmax=idesc;
    iout1=id1;
    plot(ioutmax,youtmax,'k*','markersize',10)
else
    youtmax=yasc;
    ioutmax=iasc;
    iout1=ia1;
    plot(N-ioutmax+1,youtmax,'k*','markersize',10)
end

threshvec=zeros(N,1);
threshvec(iout1(1:ioutmax))=1;
threshvec(setdiff(1:N,iout1(1:ioutmax)))=-1;



