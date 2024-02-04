function hhR=nonlinzone(x,Bright,right)
hhR=ones(length(x),1);

 DeltaR=(right-Bright);

 hhR=max(sign(Bright-x),min(1+sign(Bright+DeltaR-x), (sin(pi*(x-Bright+DeltaR)/(2*DeltaR))).^2 ));
 hhR=1-hhR;
