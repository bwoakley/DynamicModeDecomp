% Do the stats
cut=40;
t=0:.05:19.95;
Htot=zeros(2,1);
Etot=Htot; Itot=Htot;

n_str{1}=strcat(' levy');
n_str{2}=strcat(' gauss');
n_str{3}=strcat(' det');

t_str{1}=strcat('ow');
t_str{2}=strcat('ftle');
t_str{3}=strcat('meso');
t_str{4}=strcat('geo');
t_str{5}=strcat('fte');

for ni=3
    for ti=1:5
        Htot=zeros(2,1);
Etot=Htot; Itot=Htot;
        for i=1:200
            loadfile=strcat(n_str{ni},'_data/real_',t_str{ti},'_',num2str(i),'.mat')
            sss=strcat('load',loadfile);
            eval(sss);
        
            for j=1:32768
                P=xz(j,1:cut); P=P(:);
                Pdiff=diff(P);
                ind=find(Pdiff~=0);
                if(length(ind)~=0)
                    if(P(1)>.1)
                        Htot(2)=Htot(2)+1;
                    elseif(P(1)<-.1)
                        Etot(2)=Etot(2)+1;
                    else
                        Itot(2)=Itot(2)+1;
                    end
                else
                    if(P(end)>.1)
                        Htot(1)=Htot(1)+1;
                    elseif(P(end)<-.1)
                        Etot(1)=Etot(1)+1;
                    else
                        Itot(1)=Itot(1)+1;
                    end
                end
            end
        end
        
        Hr=Htot(1)/sum(Htot);
        Er=Etot(1)/sum(Etot);
        Ir=Itot(1)/sum(Itot);

        savefile=strcat(n_str{ni},'_',t_str{ti},'table_results.mat');
        sss=strcat('save',savefile,' Hr Er Ir');
        eval(sss); 
    end
end