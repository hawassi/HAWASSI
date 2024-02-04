% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function u = funOprt_UgExact(k,d,g,~)
J=length(d);
u = zeros(length(k),J);
for j=1:J
    for i=1:length(k)
        if k(i)==0
            u(i,j)=sqrt(g*d(j));
        else
            u(i,j) = sign(k(i))*sqrt(g)/2./(k(i).*tanh(d(j).*k(i)))^(1/2).*(tanh(d(j)*k(i))+k(i)*(1-tanh(d(j)*k(i))^2).*d(j));
        end
    end
end
end