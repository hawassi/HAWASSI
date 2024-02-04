function w = OmExact(k,d,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sign(k).*sqrt(9.81*k.*tanh(d(j)*k)); 
end;
end