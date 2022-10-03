function [W1,W2]=noisetype(KK,n_type)

if(n_type=='N')
    W1=zeros(KK,1);
    W2=W1;
elseif(n_type=='G')
    W1=sqrt(2*.000001)*randn(KK,1);
    W2=sqrt(2*.0000001)*randn(KK,1);
else
    X=stblrnd(.75,0,1,0,KK,1)';
    Phi=unifrnd(0,2*pi,KK,1)';
    W1=X.*cos(Phi); W1=.0001*W1(:);
    W2=X.*sin(Phi); W2=.0001*W2(:);
end