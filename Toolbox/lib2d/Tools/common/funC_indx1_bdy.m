function Indx=funC_indx1_bdy(MM)
%[GMx,GMy]=gradient(MM);

[GMFx,GMFy]=funOprt_FDGradient(MM,1,1,1);
[GMBx,GMBy]=funOprt_FDGradient(MM,1,1,-1);
GMFx(GMFx>0)=0;GMFy(GMFy>0)=0;
GMBx(GMBx<0)=0;GMBy(GMBy<0)=0;
GMx=(GMFx+GMBx);
GMy=(GMFy+GMBy);

bb=(abs(GMx)+abs(GMy));%.*(1-MM);
Indx=find(bb>0);
end