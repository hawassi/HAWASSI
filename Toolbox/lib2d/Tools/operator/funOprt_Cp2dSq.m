function Cp2dSqA_hat=funOprt_Cp2dSq(Cp2dSqInterp,Z_hat,Z)
Cp2dSqmin=Cp2dSqInterp.Opmin;
Cp2dSqplus=Cp2dSqInterp.Opplus;
gammin  =Cp2dSqInterp.gammin;
gamplus =Cp2dSqInterp.gamplus;

if isfield(Cp2dSqInterp,'gammid')
gammid=Cp2dSqInterp.gammid;   
Cp2dSqmid=Cp2dSqInterp.Opmid;
Cp2dSqA_hat=(funOprt_HSS2d(Cp2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dSqplus,gamplus,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dSqmid,gammid,Z_hat,Z));    
else
Cp2dSqA_hat=(funOprt_HSS2d(Cp2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dSqplus,gamplus,Z_hat,Z)); 
end