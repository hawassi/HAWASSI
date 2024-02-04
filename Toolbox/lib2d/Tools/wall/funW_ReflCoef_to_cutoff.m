function CutOff = funW_ReflCoef_to_cutoff(R,nonlin)
% %Poly=[2.017 -0.4951 0.1819 -0.0001131-R];
% Poly=[2.017 -0.4951 0.1819 -R];
% root=abs(roots(Poly));
% CutOff=root(1);
if nonlin==1
a=0.8086;%0.6733;
b=0.2413;%0.3597;
c=0.0044;%0.00695;
if R==1
CutOff=1;
elseif R==0
CutOff=0;   
else
CutOff=(-b+sqrt(b^2-4*a*(c-R)))/(2*a); 
end
else
a=0.7289;%0.6752;
b=0.3229;%0.4766;
c=-0.00058;%0.005786;
if R==0
CutOff=0;
else
CutOff=(-b+sqrt(b^2-4*a*(c-R)))/(2*a);  
end
end
