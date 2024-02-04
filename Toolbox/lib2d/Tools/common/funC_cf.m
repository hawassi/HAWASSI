 function cfun=funC_cf(xi,Xi,del)
cfun = funC_Heaviside(xi-Xi).*max(sign(xi-(Xi+del)),(1-cos((xi-Xi)*pi/del))/2);
