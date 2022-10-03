function [bu1,bs,bv]=rsvd(Ai,q,p)
n=size(Ai,2);
kk=floor(n/p);
if(q~=0)
    [Q,~]=qr((Ai*Ai').^q*Ai*randn(n,kk));
else
    [Q,~]=qr(Ai*randn(n,kk));
end
Q=Q(:,1:kk);

% ind_1=find(isnan(Ai(:)) || isinf(Ai(:)));
% Ai(ind_1)=0;
save tmp2.mat Ai Q
[buhat,bs,bv]=svd(Q'*Ai);
bu1=Q*buhat;