function L2dA_hat=funOprt_L2d_runup(g,Om2dsqInterp,Z_hat,Z)

L2dA_hat=(funOprt_HSS2d(Om2dsqInterp.Opmin1,Om2dsqInterp.gammin1,Z_hat,Z)+...
         funOprt_HSS2d(Om2dsqInterp.Opplus1,Om2dsqInterp.gamplus1,Z_hat,Z)+...
         funOprt_HSS2d(Om2dsqInterp.Opmid1,Om2dsqInterp.gammid1,Z_hat,Z)+....
         funOprt_HSS2d(Om2dsqInterp.Opmin2,Om2dsqInterp.gammin2,Z_hat,Z)+...
         funOprt_HSS2d(Om2dsqInterp.Opplus2,Om2dsqInterp.gamplus2,Z_hat,Z)+...
         funOprt_HSS2d(Om2dsqInterp.Opmid2,Om2dsqInterp.gammid2,Z_hat,Z))./g;    
