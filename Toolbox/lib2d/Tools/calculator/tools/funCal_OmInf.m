function w = OmInf(k,d) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sign(k).*sqrt(9.81*abs(k)); 
end;
end