function I = flood2(I,c,x,y)

LastFlood = zeros(size(I));
Flood = LastFlood;
Flood(x,y) = 1;
Mask = (abs(I-I(x,y))<=.1);
FloodFilter = [0,1,0; 1,1,1; 0,1,0];
while any(LastFlood(:) ~= Flood(:))
LastFlood = Flood;
Flood = conv2(Flood,FloodFilter,'same') & Mask;
Flood=double(Flood);
end
I(find(Flood)) = c;