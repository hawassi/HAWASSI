 function cfun=cf(xi,Xi,del)
% xi=-7:0.1:9;
% Xi=0;
% del=4;
cfun = Heaviside(xi-Xi).*max(sign(xi-(Xi+del)),(1-cos((xi-Xi)*pi/del))/2);
% cf1 = cf.*max(sign(xi-(Xi+del)),(sin((xi-Xi)*pi/2/del).^2));
% % cf = cf.*
% figure(2)
% plot(xi,cf)