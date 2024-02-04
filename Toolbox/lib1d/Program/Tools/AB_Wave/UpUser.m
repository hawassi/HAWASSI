function u = UpUser(k,d,funDisp)
J=length(d);
u = zeros(length(k),J);
O=OmUser(k,d,funDisp);
for j=1:J
    for i=1:length(k)
        if k(i)==0
            u(i,j)=sqrt(9.81*d(j));
        else
            u(i,j) = O(i,j)./k(i);
        end
    end
end
end