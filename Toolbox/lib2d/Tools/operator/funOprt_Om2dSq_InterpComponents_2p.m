function Om2d_sq=funOprt_Om2dSq_InterpComponents_2p(g,Oprt_Om2d_fun,InterpDisp,dom,current)
gammin=InterpDisp.gam_min;
gamplus=InterpDisp.gam_plus;
D_min=InterpDisp.D_min;
D_plus=InterpDisp.D_plus;
if current.check==0
Om2d_sq.Opmin=Oprt_Om2d_fun(dom.Kx,dom.Ky,D_min,g,[1 1]).^2;
Om2d_sq.Opplus=Oprt_Om2d_fun(dom.Kx,dom.Ky,D_plus,g,[1 1]).^2;
else
Om2d_sq.Opmin=Oprt_Om2d_fun(dom.Kx,dom.Ky,D_min,g,[1 1],current.ux,current.uy).^2;
Om2d_sq.Opplus=Oprt_Om2d_fun(dom.Kx,dom.Ky,D_plus,g,[1 1],current.ux,current.uy).^2;    
end
Om2d_sq.gammin=gammin;
Om2d_sq.gamplus=gamplus;

