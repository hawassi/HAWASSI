% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function u = funOprt_UgSWE(k,d,g,~)
J=length(d);
u = zeros(length(k),J);
for j=1:J
            u(:,j)=sqrt(g*d(j));
end
end