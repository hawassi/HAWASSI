function w = funOprt_OmBBM(k,d,g,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sqrt(g*d(j))*k./(1+(d(j)*k).^2/6); 
end;
end