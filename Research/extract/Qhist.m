% load Qhistory.mat;
tic
t=.2:.2:20;
Etot=zeros(length(t),1);
Itot=Etot;
Htot=Etot;
for ri=1:Ri
    ri
    Qhis(:,:,:)=Qa(:,:,ri,:);
ind=find(Qhis>1);Qhis(ind)=1;
ind=find(Qhis<-1);Qhis(ind)=-1;
ind=find((Qhis<1)&(Qhis>-1));Qhis(ind)=0;
for i=1:k
    for j=1:k
        QQQ=Qhis(i,j,:);QQQ=QQQ(:);
        Qdiff=diff(QQQ);
        ind=find(Qdiff~=0);
        if(length(ind)~=0)
        lengthy=diff(ind);
        lengthy=[ind(1); lengthy];
        for kk=1:length(ind)
           if(QQQ(ind(kk))>0.1)
               Htot(lengthy(kk))=Htot(lengthy(kk))+1;
           elseif(QQQ(ind(kk))<-0.1)
               Etot(lengthy(kk))=Etot(lengthy(kk))+1;
           else
               Itot(lengthy(kk))=Itot(lengthy(kk))+1;
           end
        end
        else
           if(QQQ(end)>0.1)
               Htot(end)=Htot(end)+1;
           elseif(QQQ(end)<-0.1)
               Etot(end)=Etot(end)+1;
           else
               Itot(end)=Itot(end)+1;
           end
        end
    end
end

end
Htot=Htot/sum(Htot);
Itot=Itot/sum(Itot);
Etot=Etot/sum(Etot);
subplot(1,2,1)
loglog(t,Htot)
hold on;
loglog(t,Itot,'r')
loglog(t,Etot,'k')
legend
subplot(1,2,2)
semilogy(t,Htot)
hold on;
semilogy(t,Itot,'r')
semilogy(t,Etot,'k')
legend
toc