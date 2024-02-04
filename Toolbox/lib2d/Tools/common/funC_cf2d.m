 function cfun=funC_cf2d(xi,Xi,del,N)
cfun = funC_Heaviside(xi-Xi).*max(sign(xi-(Xi+del)),(1-cos((xi-Xi)*pi/del))/2);
cfun =repmat(cfun,N,1);
