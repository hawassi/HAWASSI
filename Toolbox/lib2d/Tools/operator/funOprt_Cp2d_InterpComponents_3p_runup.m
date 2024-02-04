function Cp2d=funOprt_Cp2d_InterpComponents_3p_runup(g,bath,Oprt_Cp2d_fun,InterpDisp,dom,current)

gammin1=InterpDisp.sp_gam_min1;
gamplus1=InterpDisp.sp_gam_plus1;
gammid1=InterpDisp.sp_gam_mid1;
D_min1=InterpDisp.H_min1;
D_mid1=InterpDisp.H_mid1;
D_plus1=InterpDisp.H_plus1;

kkx = dom.Kx;
kky = dom.Ky;

if current.check==0
    Cp2d.Opmin1=Oprt_Cp2d_fun(kkx,kky,D_min1,g);
    Cp2d.Opmid1=Oprt_Cp2d_fun(kkx,kky,D_mid1,g);
    Cp2d.Opplus1=Oprt_Cp2d_fun(kkx,kky,D_plus1,g);
else
    Cp2d.Opmin1=Oprt_Cp2d_fun(kkx,kky,D_min1,g,current.ux,current.uy);
    Cp2d.Opmid1=Oprt_Cp2d_fun(kkx,kky,D_mid1,g,current.ux,current.uy);
    Cp2d.Opplus1=Oprt_Cp2d_fun(kkx,kky,D_plus1,g,current.ux,current.uy);
end

gammin2=InterpDisp.sp_gam_min2;
gamplus2=InterpDisp.sp_gam_plus2;
gammid2=InterpDisp.sp_gam_mid2;
D_min2=InterpDisp.H_min2;
D_mid2=InterpDisp.H_mid2;
D_plus2=InterpDisp.H_plus2;
if current.check==0
    Cp2d.Opmin2=Oprt_Cp2d_fun(kkx,kky,D_min2,g);
    Cp2d.Opmid2=Oprt_Cp2d_fun(kkx,kky,D_mid2,g);
    Cp2d.Opplus2=Oprt_Cp2d_fun(kkx,kky,D_plus2,g);
else
    Cp2d.Opmin2=Oprt_Cp2d_fun(kkx,kky,D_min2,g,current.ux,current.uy);
    Cp2d.Opmid2=Oprt_Cp2d_fun(kkx,kky,D_mid2,g,current.ux,current.uy);
    Cp2d.Opplus2=Oprt_Cp2d_fun(kkx,kky,D_plus2,g,current.ux,current.uy);
end

Cp2d.gammin1=gammin1;
Cp2d.gamplus1=gamplus1;
Cp2d.gammid1=gammid1;
Cp2d.gammin2=gammin2;
Cp2d.gamplus2=gamplus2;
Cp2d.gammid2=gammid2;