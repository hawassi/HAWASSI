function Int2mn=fun_IntegVarphimVarphin(m,n,D,Z)
%integral_-D_Zmin varphi_m varphi_n dz
if m==0 && n==0
    Int2mn=D+Z;
elseif m==n
    Int2mn=(D+Z)/2;
else
    Int2mn=0;
end
end