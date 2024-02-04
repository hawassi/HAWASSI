function Cp2d=funOprt_Cp2d_InterpComponents_3p(g,Oprt_Cp2d_fun,InterpDisp,dom,current)
gammin=InterpDisp.gam_min;
gamplus=InterpDisp.gam_plus;
gammid=InterpDisp.gam_mid;

D_min=InterpDisp.D_min;
D_mid=InterpDisp.D_mid;
D_plus=InterpDisp.D_plus;
if current.check==0
Cp2d.Opmin=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_min,g);
Cp2d.Opmid=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_mid,g);
Cp2d.Opplus=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_plus,g);
else
Cp2d.Opmin=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_min,g,current.ux,current.uy);
Cp2d.Opmid=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_mid,g,current.ux,current.uy);
Cp2d.Opplus=Oprt_Cp2d_fun(dom.Kx,dom.Ky,D_plus,g,current.ux,current.uy);
end
Cp2d.gammin=gammin;
Cp2d.gamplus=gamplus;
Cp2d.gammid=gammid;