%function CdC/deta (exact) 
%is used for runup model 
function detaC2 =detaC2Exact(k,d)
J=length(d);g=9.81;
detaC2 = zeros(length(k),J);
for j=1:J
            detaC2(:,j) = (g).*(1-tanh(k.*d(j)).^2);
end

                 
    