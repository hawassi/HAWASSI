function w = OmBBM(k,d,~) 
w = zeros(length(k),length(d));
for j=1:length(d)
    w(:,j) = sqrt(9.81*d(j))*k./(1+(d(j)*k).^2/6); 
end;
end