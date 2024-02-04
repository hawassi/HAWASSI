% Version May 2011
% Ug (group velocity)
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function u = funOprt_UgBBM(k,d,g,~)
J=length(d);
u = zeros(length(k),J);
for j=1:J
    for i=1:length(k)
        if k(i)==0
            u(i,j)=sqrt(g*d(j));
        else
            u(i,j) = sqrt(g*d(j))*(1-(d(j)*k(i))^2/6)/(1+(d(j)*k(i))^2/6)^2;
        end
    end
end
end