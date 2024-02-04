function [Output]=funSA_diffraction_or_radiation(par,ship,wave,bottom,calc)
Nk=par.Nk;
k0=par.k0;
D=bottom.depth;
w0=par.w0;
Dbdy=par.Dbdy;
g=par.g;
Evmodes=calc.evmodes;
shipform=par.shipform;
shipform_aprox=par.shipform_aprox;
Dxx=par.Dxx;
dxx=par.dxx;
xx=par.xx;
indxBdy=par.indxBdy;
Ndom=par.Ndom;
dxDbdy=par.dxDbdy;
dr=par.dr;
etaWl=[0 0];
dxetaWl=[0 0];
Xbdy=par.Xbdy;
chiS=par.chiS;
Xc  =par.Xc;%ship.cog(1);
Zc  =ship.cog(2);
XshipL=par.XshipL;
XshipR=par.XshipR;
dxZshipL=par.dxZshipL;
dxZshipR=par.dxZshipR;
ZshipL=par.ZshipL;
ZshipR=par.ZshipR;
AreaO=par.AreaO;
B=ship.width;

diff_Forces=zeros(Nk,7)*nan;
diff_RefTransCoef=zeros(Nk,3)*nan;

indx1=par.indx1;
indx2=par.indx2;
chiSwl=par.chiSwl;
GammaH=zeros(size(xx));
GammaP=zeros(size(xx));
dzX_L=par.dzX_L;
dzX_R=par.dzX_R;

rad_addedMass=zeros(Nk,10);
rad_dampCoef=zeros(Nk,10);
rad_waveAmpl=zeros(Nk,7);

% figure;
global Idstop


for ii=1:Nk
    IncAmpl=par.Ainc(ii);
   % disp(['lambda per D= ', num2str(2*pi./k0(ii)/D)])
    sig0=w0(ii)^2./g;
    
    k0L=invOmExact(w0(ii),Dbdy(1));
    k0R=invOmExact(w0(ii),Dbdy(end));
    
    try
        kappaL = [1i.*k0L invOmEvanescent(sig0*Dbdy(1),Evmodes)'./Dbdy(1)];
        kappaR = [1i.*k0R invOmEvanescent(sig0*Dbdy(end),Evmodes)'./Dbdy(end)];
    catch
        disp('here')
        continue;
    end
    
    DZeta=Dxx+shipform;
    
    dxDZeta=funOprt_FDGradient1d(DZeta,dxx,2);
    %         try
    %         [Oprt]=fun_operator_setup(kk,fun_Up,fun_Om,nu_p,DZeta,ProjDir,ProjSavename,modeldyn);
    %         catch
    %             continue;
    %         end
    
    %         tauG=fun_tangential_flow_global(kk,OpF0,cfSA,shipform,Nonlin);
    tauGInterp=xx;%fun_tangential_flow_global_Interp(g,kk,Oprt,shipform,cfSA,Nonlin);
    dxtauGInterp=gradient(tauGInterp,dxx);%Ifft(1i.*kk.*fft(tauGInterp));
    
    %     Fluxtaug  =ifft(funOprt_Flux(g,kk,Cp,shipform,tauG,Nonlin),'symmetric');
    %     FluxtaugInterp  =ifft(funOprt_Flux_interp(g,kk,Oprt,shipform,tauG,Nonlin),'symmetric');
    
    %     figure;
    %     subplot(2,1,1)
    %     plot(xx,tauGInterp,'b',xx,tauG,'r')
    %     subplot(2,1,2)
    %     plot(xx,FluxtaugInterp,'b',xx,Fluxtaug,'r')
    
    %     psiA=fun_radiation_potential_anzats(kk,xx,OpL0,OpF0,tauG,shipform,Idmotion,Nonlin,ChiAdj,xShip,chiSwl,1);
    %     psiAInterp=fun_radiation_potential_anzats_Interp(g,kk,xx,Oprt,tauGInterp,shipform,Idmotion,Nonlin,ChiAdj,xShip,chiSwl,1,Idplot);
    %     dxpsiAInterp=Ifft(1i.*kk.*fft(psiAInterp));
    %     GammaAB=(((D+shipform).^2+psiAInterp)./2./(D+shipform)).*chiSwl;
    %     dxGammaAB=gradient(GammaAB,dxx);
    
    %     figure;
    %     subplot(2,1,1)
    %     plot(xx(indx1:indx2),Gamma,'b',xx(indx1:indx2),GammaAB(indx1:indx2),'r')
    %     subplot(2,1,2)
    %     plot(xx(indx1:indx2),dxGamma,'b',xx(indx1:indx2),dxGammaAB(indx1:indx2),'r')
    kappaS=[kappaL.' kappaR.'];
    %     psibdy=[psiAInterp(indx1) psiAInterp(indc) psiAInterp(indx2)];
    %     dxpsibdy=[dxpsiAInterp(indx1) dxpsiAInterp(indc) dxpsiAInterp(indx2)];
    tauGbdy=tauGInterp(indxBdy);%[tauGInterp(indx1) tauGInterp(indxlc) tauGInterp(indc) tauGInterp(indxcr) tauGInterp(indx2)];
    dxtauGbdy=dxtauGInterp(indxBdy);%[dxtauGInterp(indx1) dxtauGInterp(indxlc) dxtauGInterp(indc) dxtauGInterp(indxcr) dxtauGInterp(indx2)];
    DZetabdy=DZeta(indxBdy);%[DZeta(indx1) DZeta(indxlc) DZeta(indc) DZeta(indxcr) DZeta(indx2)];
    
    dxDZetabdy=zeros(size(DZetabdy));%%dxDZeta(indxBdy);%[dxDbdy(1) dxDZeta(indxlc) dxDZeta(indc) dxDZeta(indxcr) dxDbdy(5)];
    if calc.Id==1
        funSA_diffraction_calculation;
    elseif calc.Id==2
        funSA_radiation_calculation;
    end   
    
    if mod(ii,floor(0.1*Nk))==0
    set(par.jProgressBar,'Maximum',Nk, 'Value',ii);
    par.jProgressBar.setStringPainted( true );
    ETA=remain_time(ii,Nk);
    par.statusbarObj.setText(['time=', num2str(ETA)]);
    end
    
    fun_Stop_iter_button(par.jbStop);
   
    if Idstop==1
        set(par.jProgressBar,'Maximum',Nk, 'Value',Nk);
        par.jProgressBar.setStringPainted( true );
        break;
    end
end
if Idstop==1
par.statusbarObj.setText(['Terminated.']);   
else
par.statusbarObj.setText(['']);
end

if calc.Id==1
     Output.diff.Forces=diff_Forces;
     Output.diff.ReflTrans=diff_RefTransCoef;
     
elseif calc.Id==2
    Output.rad.addedmass=rad_addedMass;
    Output.rad.dampCoef=rad_dampCoef;
    Output.rad.waveAmpl=rad_waveAmpl;
end
end