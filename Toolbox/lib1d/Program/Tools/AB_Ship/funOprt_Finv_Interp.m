function FinvVar_hat=funOprt_Finv_Interp(g,kk,Oprt,var,bathtype)
var_hat=fft(var);

if strcmpi(bathtype,'F') || strcmpi(bathtype,'B') || strcmpi(bathtype,'U')
    OneperC2min=1./Oprt.Csqmin;
    OneperC2plus=1./Oprt.Csqplus;
    gm          = Oprt.gammin;
    gp          = Oprt.gamplus;
    if Oprt.interp      ==2
        OneperC2Var_hat   = HSS(OneperC2plus,gp,var_hat,var)+HSS(OneperC2min,gm,var_hat,var);
    else
        OneperC2mid=1./Oprt.Csqmid;
        gc          = Oprt.gammid;
        OneperC2Var_hat   = HSS(OneperC2plus,gp,var_hat,var)+HSS(OneperC2mid,gc,var_hat,var)+HSS(OneperC2min,gm,var_hat,var);
    end
else
    gm1=Oprt.interp.m1;gc1=Oprt.interp.c1;gp1=Oprt.interp.p1;
    gm2=Oprt.interp.m2;gc2=Oprt.interp.c2;gp2=Oprt.interp.p2;
    
    OneperC2m1=1./Oprt.C2m1;OneperC2c1=1./Oprt.C2c1;OneperC2p1=1./Oprt.C2p1;
    OneperC2m2=1./Oprt.C2m2;OneperC2c2=1./Oprt.C2c2;OneperC2p2=1./Oprt.C2p2;
    
    OneperC2Var_hat  = Op_C2r(OneperC2m1,gm1,OneperC2p1,gp1,OneperC2c1,gc1,OneperC2m2,gm2,OneperC2p2,gp2,OneperC2c2,gc2,var,var_hat);
end
FinvVar_hat=g*OneperC2Var_hat./(1i.*kk);
FinvVar_hat(1)=0;