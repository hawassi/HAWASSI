function CutOff = ReflCoef_to_cutoff(R)
%Poly=[2.017 -0.4951 0.1819 -0.0001131-R];
Poly=[2.017 -0.4951 0.1819 -R];
root=abs(roots(Poly));
CutOff=root(1);
end