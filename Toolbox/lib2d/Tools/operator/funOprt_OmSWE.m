function w = funOprt_OmSWE(k,d,g,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sqrt(g*d(j)).*k; 
end;
end