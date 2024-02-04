%Bore to hydraulic jump condition
function Op_var_hat = Jump_Cond(u,eta,depth,Op_var,var,flag,k,x,Chi)
H=depth+eta.*Chi;
g=9.81;B=0.5;
dxHu=gradient(H.*u,x(2)-x(1));

Cond1=((abs(dxHu)+dxHu)/2)./sqrt(g.*H);
Cond2=B;

% Cond1=H./depth;
% Cond2=0.82;


if flag==0
Op_var(Cond1>Cond2)=sqrt(g.*H(Cond1>Cond2)).*var(Cond1>Cond2);%C_var
else
Op_var(Cond1>Cond2)=(g./2)./sqrt(g.*H(Cond1>Cond2)).*var(Cond1>Cond2);%Cder_var    
end
Op_var_hat=fft(Op_var);
end