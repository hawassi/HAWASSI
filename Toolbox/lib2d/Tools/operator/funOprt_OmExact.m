function w = funOprt_OmExact(k,d,g,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sign(k).*sqrt(g*k.*tanh(d(j)*k)); 
end
end