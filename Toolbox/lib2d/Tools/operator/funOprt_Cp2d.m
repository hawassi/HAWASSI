function Cp2dA_hat=funOprt_Cp2d(Cp2dInterp,Z_hat,Z)
Cp2dmin=Cp2dInterp.Opmin;
Cp2dplus=Cp2dInterp.Opplus;
gammin  =Cp2dInterp.gammin;
gamplus =Cp2dInterp.gamplus;

if isfield(Cp2dInterp,'gammid')
gammid=Cp2dInterp.gammid;   
Cp2dmid=Cp2dInterp.Opmid;
Cp2dA_hat=(funOprt_HSS2d(Cp2dmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dplus,gamplus,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dmid,gammid,Z_hat,Z));    
else
Cp2dA_hat=(funOprt_HSS2d(Cp2dmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dplus,gamplus,Z_hat,Z)); 
end