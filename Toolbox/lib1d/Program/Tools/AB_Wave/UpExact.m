function u = UpExact(k,d,~)
J=length(d);
u = zeros(length(k),J);
O=OmExact(k,d);
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

% O=OmExact(k,d);
% Nd=length(d);Nk=length(k);
% k=repmat(k,1,Nd);
% d=repmat(d.',Nk,1);
% u=O./k;
% u(k==0)=sqrt(9.81*d(k==0));