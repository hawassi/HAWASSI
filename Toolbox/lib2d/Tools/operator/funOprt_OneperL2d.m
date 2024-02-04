function OneperL2dA_hat=funOprt_OneperL2d(g,OneperOm2dSqInterp,Z_hat,Z)
OneperOm2dSqmin=OneperOm2dSqInterp.Opmin;
OneperOm2dSqplus=OneperOm2dSqInterp.Opplus;
gammin  =OneperOm2dSqInterp.gammin;
gamplus =OneperOm2dSqInterp.gamplus;

if isfield(OneperOm2dSqInterp,'gammid')
gammid=OneperOm2dSqInterp.gammid;   
OneperOm2dSqmid=OneperOm2dSqInterp.Opmid;
OneperL2dA_hat=g*(funOprt_HSS2d(OneperOm2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(OneperOm2dSqplus,gamplus,Z_hat,Z)+...
         funOprt_HSS2d(OneperOm2dSqmid,gammid,Z_hat,Z));    
else
OneperL2dA_hat=g*(funOprt_HSS2d(OneperOm2dSqmin,gammin,Z_hat,Z)+...
         funOprt_HSS2d(OneperOm2dSqplus,gamplus,Z_hat,Z)); 
end