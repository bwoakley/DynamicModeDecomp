profile on -timer real
fp_heir
profile viewer
break
%%
profile off

%%
load ../data_0/final_positions5.mat
xa=xb{1};
ya=yb{1};
for kk=1:10
    figure(1);clf;
for i=1:25
    plot(xa(:,i),ya(:,i),'.k')
    drawnow;
end
end