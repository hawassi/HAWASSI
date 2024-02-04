function OneperOm2d_sq=funOprt_OneperOm2dSq_InterpComponents_2p(g,Oprt_Om2d_fun,InterpDisp,dom,current)
gammin=InterpDisp.gam_min;
gamplus=InterpDisp.gam_plus;
D_min=InterpDisp.D_min;
D_plus=InterpDisp.D_plus;
if current.check==0
OneperOm2d_sq.Opmin=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_min,g,[1 1]).^2);
OneperOm2d_sq.Opplus=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_plus,g,[1 1]).^2);
else
OneperOm2d_sq.Opmin=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_min,g,[1 1],current.ux,current.uy).^2);
OneperOm2d_sq.Opplus=1./(Oprt_Om2d_fun(dom.Kx,dom.Ky,D_plus,g,[1 1],current.ux,current.uy).^2);    
end
OneperOm2d_sq.gammin=gammin;
OneperOm2d_sq.gamplus=gamplus;

