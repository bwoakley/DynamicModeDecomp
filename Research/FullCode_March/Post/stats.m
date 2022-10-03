% Do the stats
t_0=0;
t_f=20;
part_space=.05;

t=t_0:part_space:t_f-part_space;
buf_val=40;
field_type='FTLE';
noise_type='det';

s_dir=strcat(' ../somewhere/else/');
r_dir=strcat(' ../evenmore/elsewhere/');

%Need to change the list values to correspond to those that are used in the
%partitions.

saver=strcat('save ',s_dir,noise_type,'/',field_type,'_results.mat Tot t');
if(strcmp(field_type,'PS'))
    list=[-1,0,1,2];
elseif(strcmp(field_type,'OW'))
    list=[-1,0,1];    
elseif(strcmp(field_type,'FTLE'))
    list=[-1,0,1];    
elseif(strcmp(field_type,'MESO'))
    list=[-1,0,1];    
elseif(strcmp(field_type,'GEO'))
    list=[0,1];   
elseif(strcmp(field_type,'FTE'))
    list=[-1,0,1];
elseif(strcmp(field_type,'AI'))
    list=[-1, 0, 1, 2];
elseif(strcmp(field_type,'FTCS'))
    list=[-1, 0, 1, 2];    
else
    list=[];%What else would you possibly want??
end

for i=1:length(list)
    Tot{i}=zeros(length(t),1);
end
for i=1:Ny
    load_data=strcat('load ',r_dir,noise_type,'/',field_type,'_zone',num2str(i),'.mat')
    eval(load_data);
    for j=1:size(xz,1)
        P=xz(j,:); P=P(:);
        % The next 16 lines remove the buffer region
        if(P(1)==buf_val)
            if(indm(1)<indp(1))
                P(1:indm(1))=P(indm(1)+1);
            else
                P(1:indp(1))=P(indp(1)+1);
            end
        end
        Pdiff=diff(P);
        indp=find(Pdiff>10); flag=1;
        while(~isempty(indp) && flag<50)
            for jj=1:length(indp)
                P(indp(jj)+1)=P(indp(jj));
            end
            Pdiff=diff(P);
            indp=find(Pdiff>10);
            flag=flag+1;
        end            
        
        ind=find(Pdiff~=0);
        if(length(ind)~=0)
            lengthy=diff(ind);
            lengthy=[ind(1); lengthy];
            for kk=1:length(ind)
                ind_list=find(P(ind(kk))==list);
                Tot{ind_list}(lengthy(kk))=Tot{ind_list}(lengthy(kk))+1;
            end
        else
            ind_list=find(P(end)==list);
            Tot{ind_list}(end)=Tot{ind_list}(end)+1;
        end
    end
end
eval(saver);


%% Plotting
for i=1:length(list)
Tot{i}=Tot{i}/sum(Tot{i});
end

subplot(1,2,1)
loglog(t,Tot{1})
hold on;
loglog(t,Tot{2},'r')
loglog(t,Tot{3},'k')
legend
subplot(1,2,2)
semilogy(t,Tot{1})
hold on;
semilogy(t,Tot{2},'r')
semilogy(t,Tot{3},'k')
legend