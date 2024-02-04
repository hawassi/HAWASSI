function Int1mn=fun_diff_Integ1mn(n,kappa_m,depth,draft)
%integ from -T to 0 of Q_m Testfun_n dz, Q_m=cos(k_m (D+z))/cos(k_m D),
%Testfun_n=cos((n)*pi*(T+z)/T);

A=(kappa_m*depth+n*pi);B=(kappa_m+n*pi/draft);
C=(kappa_m*depth-n*pi);D=(kappa_m-n*pi/draft);

Int1mn=(sin(A)/B+sin(C)/D-sin(A-B*draft)/B-sin(C-D*draft)/D)/2/cos(kappa_m*depth);
end
