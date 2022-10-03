function sol = HessGaussTransform(A)
[rows,cols]=size(A);
for i=1:cols-1;
    m = eye(cols);
    [MaxVal,MaxIndex] = max(abs(A(i+1:cols,1)));
    if MaxIndex > i + 1;
        t=m(i+1,:);
        m(i+1,:)=m(MaxIndex,:);
        m(MaxIndex,:)=t;
        A = m*A*m';
    end
    m=eye(cols);
    m(i+2:cols,i+1) = -A(i+2:cols,i)/(A(i+1,i));
    mi=m;
    mi(i+2:cols,i+1)=-m(i+2:cols,i+1);
    A = m*A*mi;
    %mesh(abs(A));
end
sol = A;


    
        