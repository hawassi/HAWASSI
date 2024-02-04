function Int2mn=fun_diff_Integ2mn(m,n,D,T)
%integral_z1_z2 chi_m chi_n dz
if m==0 && n==0
    Int2mn=D-T;
elseif m==n
    Int2mn=(D-T)/2;
else
    Int2mn=0;
end
end