function [U,v_his]=picard_its(x,c,n,m,p)
U=zeros(c,n);
D=zeros(c,n);
v_his=[];
vnew=linspace(-1,1,c);
%vnew=randn(1,c);
v_his=[v_his;vnew];


for i=1:c
    for j=1:n
        D(i,j)=p(j).^m*(norm(x(:,j)-vnew(i))).^2;
    end
end

for j=1:n
    for i=1:c
        if(D(i,j)>0)
            U(i,j)=((sum((D(i,j)./(D(:,j))).^(1./(m-1))))).^(-1);
        else            
            ind=find(abs(D(:,j))<=10^-8);
            if(~isempty(ind))
                U(:,j)=0;
                U(ind,j)=1/numel(ind);
                break
            end
        end
    end
end  

%%
for k=1:n
    U(:,k)=U(:,k)/sum(U(:,k));
end

err=1;
flag=1;
while(err>10^-4 && flag<1000)
    for i=1:c
        temp=0;
        for j=1:n
            temp=temp+sum((p(j).^m*U(i,j).^m*x(:,j)));
        end
        vnew(i)=temp./sum(p'.^m.*U(i,:).^m);
    end

    if(flag~=1)
        err=(norm(vold-vnew));
    end
for i=1:c
    for j=1:n
        D(i,j)=p(j).^m*(norm(x(:,j)-vnew(i))).^2;
    end
end

for j=1:n
    for i=1:c
        if(D(i,j)>0)
            U(i,j)=((sum((D(i,j)./(D(:,j))).^(1./(m-1))))).^(-1);
        else            
            ind=find(abs(D(:,j))<=10^-8);
            if(~isempty(ind))
                U(:,j)=0;
                U(i,j)=1;
                break
            end
        end
    end
end  
for k=1:n
    U(:,k)=U(:,k)/sum(U(:,k));
end

flag=flag+1;
v_his=[v_his;vnew];
vold=vnew;
end