% Interpolation C based on 3 points different depth.
function C = Op_C(Cm,gm,Cp,gp,Cc,gc,iftA,A)
    C=(HSS(Cm,gm,A,iftA)+HSS(Cp,gp,A,iftA)+HSS(Cc,gc,A,iftA));
end