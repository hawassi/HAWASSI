function Linv = OpLinv(k,d)
J=length(d);
Linv = zeros(length(k),J);
O=OmExact(k,d);
for j=1:J
    for i=1:length(k)
        if k(i)==0
            Linv(i,j)=1/d(j);
        else
            Linv(i,j) = 9.811./(O(i,j).^2);
        end
    end
end
end