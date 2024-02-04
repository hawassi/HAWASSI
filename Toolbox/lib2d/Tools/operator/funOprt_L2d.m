function L2dA_hat=funOprt_L2d(g,Om2dSqInterp,Z_hat,Z)
Om2dSqmin=Om2dSqInterp.Opmin;
Om2dSqplus=Om2dSqInterp.Opplus;
gammin  =Om2dSqInterp.gammin;
gamplus =Om2dSqInterp.gamplus;

if isfield(Om2dSqInterp,'gammid')
gammid=Om2dSqInterp.gammid;   
Om2dSqmid=Om2dSqInterp.Opmid;
L2dA_hat=(funOprt_HSS2d(Om2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Om2dSqplus,gamplus,Z_hat,Z)+...
         funOprt_HSS2d(Om2dSqmid,gammid,Z_hat,Z))/g;    
else
L2dA_hat=(funOprt_HSS2d(Om2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Om2dSqplus,gamplus,Z_hat,Z))/g; 
end