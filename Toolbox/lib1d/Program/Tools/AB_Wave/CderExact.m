%function dC/deta (exact) 
%is used for runup model 
function Cprime =CderExact(k,d)
J=length(d);g=9.81;
Cprime = zeros(length(k),J);
for j=1:J
            Cprime(:,j) = (g./2).* sqrt(k./(g.*tanh(k.*d(j)))).*(1-tanh(k.*d(j)).^2);
            if k(1)==0
            Cprime(1,j) = (g./2)./sqrt(g.*d(j));
            end
end

                 
    