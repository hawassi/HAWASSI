function Indx=funC_indx_bdy(MM)
[GMx,GMy]=gradient(MM);
bb=(abs(GMx)+abs(GMy)).*(1-MM);
Indx=find(bb>0);
end