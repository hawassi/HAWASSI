% Interpolation C based on 3 points different depth. for runup case
function C2 = Op_C2r(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,iftA,A)
    C2=(HSS(Cm1,gm1,A,iftA)+HSS(Cp1,gp1,A,iftA)+HSS(Cc1,gc1,A,iftA))+...
        (HSS(Cm2,gm2,A,iftA)+HSS(Cp2,gp2,A,iftA)+HSS(Cc2,gc2,A,iftA));
end