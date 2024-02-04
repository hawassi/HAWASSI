function OneperOm2d_sq=funOprt_OneperOm2dSq_InterpComponents_3p(g,Oprt_Om2d_fun,InterpDisp,dom)
gammin=InterpDisp.gam_min;
gamplus=InterpDisp.gam_plus;
gammid=InterpDisp.gam_mid;

D_min=InterpDisp.D_min;
D_mid=InterpDisp.D_mid;
D_plus=InterpDisp.D_plus;
OneperOm2d_sq.Opmin=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_min,g,[1 1]).^2);
OneperOm2d_sq.Opmid=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_mid,g,[1 1]).^2);
OneperOm2d_sq.Opplus=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_plus,g,[1 1]).^2);
OneperOm2d_sq.gammin=gammin;
OneperOm2d_sq.gamplus=gamplus;
OneperOm2d_sq.gammid=gammid;