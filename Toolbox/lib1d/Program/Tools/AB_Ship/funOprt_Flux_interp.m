function Fvar_hat=funOprt_Flux_interp(g,kk,Oprt,zeta,var,OrdeOp,bathtype)
dxvar_hat=1i.*kk.*fft(var);
dxvar=Ifft(dxvar_hat);

if strcmpi(bathtype,'F') || strcmpi(bathtype,'B') || strcmpi(bathtype,'U')
    
    C2min=Oprt.Csqmin;
    C2mid=Oprt.Csqmid;
    C2plus=Oprt.Csqplus;
    gm          = Oprt.gammin;
    gc          = Oprt.gammid;
    gp          = Oprt.gamplus;
    
    Fvar_hat   = (HSS(C2plus,gp,dxvar_hat,dxvar)+HSS(C2mid,gc,dxvar_hat,dxvar)+HSS(C2min,gm,dxvar_hat,dxvar))/g;
else
    gm1=Oprt.interp.m1;gc1=Oprt.interp.c1;gp1=Oprt.interp.p1;
    gm2=Oprt.interp.m2;gc2=Oprt.interp.c2;gp2=Oprt.interp.p2;
    
    C2m1=Oprt.C2m1;C2c1=Oprt.C2c1;C2p1=Oprt.C2p1;
    C2m2=Oprt.C2m2;C2c2=Oprt.C2c2;C2p2=Oprt.C2p2;
    Fvar_hat = Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,dxvar,dxvar_hat);
end

% if OrdeOp~=1
%     F1var_hat=funOp_F1(F0,zeta,kk,var);
%     Fvar_hat=Fvar_hat+ F1var_hat;
%     if OrdeOp~=2
%         F2var_hat=funOp_F2(F0,zeta,kk,var);
%         Fvar_hat=Fvar_hat+ F2var_hat;
%     end
% end