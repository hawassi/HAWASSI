function Int0mn=fun_diff_Integ0mn(m,n,kappa_m,kappa_n,D,draft)
%integ from -D to -T of Q_m Q_n dz, Q_m=cos(k_m (D+z))/cos(k_m D),

if m==n
 Int0mn=((D-draft)+sin(2*kappa_m*(D-draft))/2/kappa_m)/2/cos(kappa_m*D)/cos(kappa_m*D);
   
else

A=(kappa_m+kappa_n);
B=(kappa_m-kappa_n);

Int0mn=(sin(A*(D-draft))/A+sin(B*(D-draft))/B)/2/cos(kappa_m*D)/cos(kappa_n*D);   
end
end