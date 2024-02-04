function Lvar_hat=funOprt_L_interp(g,kk,Oprt,zeta,var,OrdeOp,bathtype)

Fvar_hat=funOprt_Flux_interp(g,kk,Oprt,zeta,var,OrdeOp,bathtype);

Lvar_hat   = -1i.*kk.*Fvar_hat;

