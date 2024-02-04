function Cp2d_sq=funOprt_Cp2dSq_InterpComponents_2p(g,Oprt_Cp2d_fun,InterpDisp,dom,current)
gammin=InterpDisp.gam_min;
gamplus=InterpDisp.gam_plus;
D_min=InterpDisp.D_min;
D_plus=InterpDisp.D_plus;
if current.check==0
Cp2d_sq.Opmin=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_min,g).^2;
Cp2d_sq.Opplus=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_plus,g).^2;
else
Cp2d_sq.Opmin=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_min,g,current.ux,current.uy).^2;
Cp2d_sq.Opplus=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_plus,g,current.ux,current.uy).^2;     
end
Cp2d_sq.gammin=gammin;
Cp2d_sq.gamplus=gamplus;

