function u = funOprt_UpSWE(k,d,g,~)
J=length(d);
u = zeros(length(k),J);
O=funOprt_OmSWE(k,d,g);
for j=1:J
    for i=1:length(k)
        if k(i)==0
            u(i,j)=sqrt(g*d(j));
        else
            u(i,j) = O(i,j)./k(i);
        end
    end
end
end