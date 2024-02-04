function Cp2dA_hat=funOprt_Cp2d_runup(Cp2dInterp,Z_hat,Z)

Cp2dmin1=Cp2dInterp.Opmin1;
Cp2dplus1=Cp2dInterp.Opplus1;  
Cp2dmid1=Cp2dInterp.Opmid1;
gammin1  =Cp2dInterp.gammin1;
gamplus1 =Cp2dInterp.gamplus1;
gammid1 =Cp2dInterp.gammid1; 

Cp2dmin2=Cp2dInterp.Opmin2;
Cp2dplus2=Cp2dInterp.Opplus2;  
Cp2dmid2=Cp2dInterp.Opmid2;
gammin2  =Cp2dInterp.gammin2;
gamplus2 =Cp2dInterp.gamplus2;
gammid2 =Cp2dInterp.gammid2; 

Cp2dA_hat=(funOprt_HSS2d(Cp2dmin1,gammin1,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dplus1,gamplus1,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dmid1,gammid1,Z_hat,Z)+....
         funOprt_HSS2d(Cp2dmin2,gammin2,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dplus2,gamplus2,Z_hat,Z)+...
         funOprt_HSS2d(Cp2dmid2,gammid2,Z_hat,Z));    
