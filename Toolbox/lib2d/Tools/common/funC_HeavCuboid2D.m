function HeavCuboid2D=funC_HeavCuboid2D(X,Y,ll,rr,bb,tt)
HeavY=funC_Heav2d(Y,bb,tt,length(X));
HeavX=funC_Heav2d(X,ll,rr,length(Y))';
HeavCuboid2D=(HeavX.*HeavY);
