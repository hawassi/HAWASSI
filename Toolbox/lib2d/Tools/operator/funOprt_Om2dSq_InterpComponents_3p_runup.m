function Om2d_sq=funOprt_Om2dSq_InterpComponents_3p_runup(g,bath,Oprt_Om2d_fun,InterpDisp,dom,current)

gammin1=InterpDisp.sp_gam_min1;
gamplus1=InterpDisp.sp_gam_plus1;
gammid1=InterpDisp.sp_gam_mid1;
D_min1=InterpDisp.H_min1;
D_mid1=InterpDisp.H_mid1;
D_plus1=InterpDisp.H_plus1;

kkx = dom.Kx;
kky = dom.Ky;

if current.check==0
    Om2d_sq.Opmin1=Oprt_Om2d_fun(kkx,kky,D_min1,g,[1 1]).^2;
    Om2d_sq.Opmid1=Oprt_Om2d_fun(kkx,kky,D_mid1,g,[1 1]).^2;
    Om2d_sq.Opplus1=Oprt_Om2d_fun(kkx,kky,D_plus1,g,[1 1]).^2;
else
    Om2d_sq.Opmin1=Oprt_Om2d_fun(kkx,kky,D_min1,g,[1 1],current.ux,current.uy).^2;
    Om2d_sq.Opmid1=Oprt_Om2d_fun(kkx,kky,D_mid1,g,[1 1],current.ux,current.uy).^2;
    Om2d_sq.Opplus1=Oprt_Om2d_fun(kkx,kky,D_plus1,g,[1 1],current.ux,current.uy).^2;
end

gammin2=InterpDisp.sp_gam_min2;
gamplus2=InterpDisp.sp_gam_plus2;
gammid2=InterpDisp.sp_gam_mid2;
D_min2=InterpDisp.H_min2;
D_mid2=InterpDisp.H_mid2;
D_plus2=InterpDisp.H_plus2;
if current.check==0
    Om2d_sq.Opmin2=Oprt_Om2d_fun(kkx,kky,D_min2,g,[1 1]).^2;
    Om2d_sq.Opmid2=Oprt_Om2d_fun(kkx,kky,D_mid2,g,[1 1]).^2;
    Om2d_sq.Opplus2=Oprt_Om2d_fun(kkx,kky,D_plus2,g,[1 1]).^2;
else
    Om2d_sq.Opmin2=Oprt_Om2d_fun(kkx,kky,D_min2,g,[1 1],current.ux,current.uy).^2;
    Om2d_sq.Opmid2=Oprt_Om2d_fun(kkx,kky,D_mid2,g,[1 1],current.ux,current.uy).^2;
    Om2d_sq.Opplus2=Oprt_Om2d_fun(kkx,kky,D_plus2,g,[1 1],current.ux,current.uy).^2;
end

Om2d_sq.gammin1=gammin1;
Om2d_sq.gamplus1=gamplus1;
Om2d_sq.gammid1=gammid1;
Om2d_sq.gammin2=gammin2;
Om2d_sq.gamplus2=gamplus2;
Om2d_sq.gammid2=gammid2;