% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
%%
function u = UgSWE(k,d,~)
J=length(d);
u = zeros(length(k),J);
for j=1:J
            u(:,j)=sqrt(9.81*d(j));
end
end