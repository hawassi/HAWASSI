function al=alias(k,frac)
n=length(k);
al=zeros(n,1);
k0=k(floor(n/frac));
al=(sign(k0-abs(k))+1)/2;
% figure;
% plot(k,al)
end