function hh=dampzone(x,Fbdy)
Bleft        =x(1)+Fbdy;
Bright       =x(end)-Fbdy;
left_bdy     =x(1);
right_bdy    =x(end);

hh=ones(1,length(x));
hhR=ones(1,length(x));
hhL=ones(1,length(x));

 DeltaR=(right_bdy-Bright);
 DeltaL=-(left_bdy-Bleft);

 hhR=max(sign(Bright-x),min(1+sign(Bright+DeltaR-x), (sin(pi*(x-Bright+DeltaR)/(2*DeltaR))).^2 ));
 hhL=min(1+sign(x-Bleft+DeltaL),max(sign(x-Bleft), (sin(pi*(x-Bleft+DeltaL)/(2*DeltaL))).^2));
 hh=hhR.*hhL;
 hh=(1-hh);



