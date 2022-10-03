


function z = denecke(x,y,p)

% This function takes as arguments two equal-length vectors x and y. p is a
% vector of length 2. The function determines if the point p lies above (+1), 
% on (0) or below (-1) the function y = f(x). 

ux = min(x(x > p(1))); lx = max(x(x < p(1)));
uy = y(x == ux); ly = y(x == lx);

if lx ~= ux
   s = spline([lx ux],[ly uy],p(2));
else
   s = uy;
end

z = sign(p(2) - s);
