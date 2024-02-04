%function dC/deta (BBM) 
%is used for runup model 
function Cprime =CderBBM(k,d)
J=length(d);g=9.81;
Cprime = zeros(length(k),J);
for j=1:J
            C0          = sqrt(g.*d(j));
            Cprime(:,j) = (g./2)./C0./(1+(k.*d(j)).^2/6)-...
                           C0.*(k.*d(j)).*k./(3.*(1+(k.*d(j)).^2/6).^2);
end