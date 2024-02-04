%function dC/deta (SWE) 
%is used for runup model 
function Cprime =CderSWE(k,d)
J=length(d);g=9.81;
Cprime = zeros(length(k),J);
for j=1:J
            Cprime(:,j) = (g./2)./sqrt(g.*d(j));
end
