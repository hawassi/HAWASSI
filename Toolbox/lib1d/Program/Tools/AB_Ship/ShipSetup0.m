%%%%%%%Ship-configuration%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%Ship-Shape%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shippar.Nship=length(shippar.data(:,1));

FlagAllSFixed=1;
for i=1:shippar.Nship
    try
        SXc(i)=str2num(cell2mat(shippar.data(i,5)));
    catch
        SXc(i)=cell2mat(shippar.data(i,5));
    end
    try
        SZc(i)=str2num(cell2mat(shippar.data(i,6)));
    catch
        SZc(i)=cell2mat(shippar.data(i,6));
    end
    try
        SThetac(i)=str2num(cell2mat(shippar.data(i,7)));
    catch
        SThetac(i)=cell2mat(shippar.data(i,7));
    end
    if  ~strcmpi(shippar.data(i,2),'Fixed')
    FlagAllSFixed=0;
    end
   
end




[shippar.form]= shapeformsetup(x,shippar,SXc,SZc,SThetac,bath.name);


nup=influx.nu_p;
gradDepth=gradient(-par.bathy,par.dx);
kappaS=zeros((shippar.Evmode+1),shippar.Nship*2);
Dwl=zeros(shippar.Nship,2);
dxDwl=zeros(shippar.Nship,2);
iterx=1;
shipform=shippar.form.shapeXcZc(end,:).';
for jj=1:shippar.Nship
    indXc=closest(par.x,shippar.form.xShip(jj,2));
    SXc=par.x(indXc);
    indXL=closest(par.x,shippar.form.xShip(jj,1));
    indXR=closest(par.x,shippar.form.xShip(jj,3));
    DL=-par.bathy(indXL);DR=-par.bathy(indXR);
    dxDL=gradDepth(indXL);dxDR=gradDepth(indXR);
    kpL=invOmExact(nup,DL);
    kpR=invOmExact(nup,DR);
    sig0=nup.^2./g;
    if shippar.Evmode>=1
        kappaDL=invOmEvanescent(sig0*DL,shippar.Evmode)';
        kappaDR=invOmEvanescent(sig0*DR,shippar.Evmode)';
        kappaL=[-1i*kpL kappaDL./DL].';
        kappaR=[-1i*kpR kappaDR./DR].';
       
        kappaS(:,iterx:iterx+1)=[kappaL kappaR];
    else
        kappaS(1,iterx:iterx+1)=[-1i*kpL -1i*kpR];
    end
    Dwl(jj,:)=[DL DR];
    dxDwl(jj,:)=[dxDL dxDR];
    iterx=iterx+2;
    
    shipform(indXL)=shipform(indXL+1); %%% % due to the barge make sure the flux at ship bottom at water line is not zero.
    shipform(indXR)=shipform(indXR-1);
end
shippar.kappa=kappaS;
shippar.Dwl=Dwl;
shippar.dxDwl=dxDwl;
assignin('base','shippar',shippar);

chiS0=par.cfSA.';



shippar.tangenflowGlobal=par.x.';%fun_tangential_flow_global(k,Oprt.F0,chiS0,shipform,model.nonlinear);%
%%%%%%%%%%%%%%%%%%%
% figure;
% subplot(2,1,1)
% plot(x,shippar.tangenflowGlobal,'r');
% subplot(2,1,2)
% Fluxtaug  =ifft(funOprt_Flux(g,k,Oprt.Upd,shipform,shippar.tangenflowGlobal,1),'symmetric');
% plot(x,Fluxtaug,'r');


A=0.1;t0=107;
eta=A*sin(influx.k_p.*x-influx.nu_p*t0).*par.cfSA;
phihat=fft(eta)./(1i.*OmExact(k,par.depth)).';
phihat(abs(phihat)==Inf)=0;
phi=g*Ifft(phihat).*par.cfSA;

chiSWl=shippar.form.chiSl(:,end)+shippar.form.chiWl(:,end);

if shippar.Nship==1
     [phif,tau,Fs,FEv,Lev,beta,dxi_Kdiff,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_Ndom...
        (shippar.calclin,k,x.',shippar.form.xShip,shippar.form.Sdraft,shippar.form.Slength,0,phi.',eta.',shippar.form.chi,shippar.form.chiWl,...
        -par.bathy,shippar.Evmode,shippar.form.XZrot0,influx.lambda_p,shippar.kappa,shippar.tangenflowGlobal,Oprt.L0,Oprt.F0,model.nonlinear,...
        influx.gen.nonlinAdj,shippar.form.shapeXcZc,shippar.form.shapedat);
%     [phif0,tau0,Fs0,FEv0,Lev0,beta0,dxi_Kdiff0,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_1...
%         (shippar.calclin,k,x.',shippar.form.xShip,shippar.form.Sdraft,shippar.form.Slength,0,phi.',eta.',shippar.form.chi,shippar.form.chiWl,...
%         -par.bathy,shippar.Evmode,shippar.form.XZrot0,influx.lambda_p,shippar.kappa,shippar.tangenflowGlobal,Oprt.L0,Oprt.F0,model.nonlinear,influx.gen.nonlinAdj);
%     
%     [phif,Fs,FEv,Lev,beta,dxi_Kdiff,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3...
%         (shippar.calclin,k,x,shippar.form.xShip,shippar.form.Sdraft,shippar.form.Slength,0,phi.',eta.',shippar.form.chi,shippar.form.chiWl,...
%         -par.bathy.',shippar.Nship,shippar.Evmode,shippar.kappa,0,0,shippar.form.nutild,shippar.form.XZrot0(1,:));
   
% figure;
% plot(x,phif0,'b',x,phif,'--r')
% betaa=[beta0.x beta0.z beta0.theta beta.x beta.z beta.theta]
% Fss=[max(Fs0) max(Fs)]
% Fevs=[max(FEv0) max(FEv)]

elseif shippar.Nship==2
   
    [phif,tau,Fs,FEv,Lev,beta,dxi_Kdiff,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_2barges...
        (shippar.calclin,k,x.',shippar.form.xShip,shippar.form.Sdraft,shippar.form.Slength,0,phi.',eta.',shippar.form.chi,shippar.form.chiWl,...
        -par.bathy,shippar.Evmode,shippar.form.XZrot0,influx.lambda_p,shippar.kappa,shippar.tangenflowGlobal,Oprt.L0,Oprt.F0,model.nonlinear,influx.gen.nonlinAdj);
end

FluxF0  =ifft(funOprt_Flux(g,k,Oprt.Upd,eta.',phi.',1),'symmetric');
FluxF  =ifft(funOprt_Flux(g,k,Oprt.Upd,eta.',phif,1),'symmetric');

Lphif0=Ifft(funOprt_L(k,Oprt.L0,eta.',phi.',1,ChiAdj));
Lphif=Ifft(funOprt_L(k,Oprt.L0,eta.',phif,1,ChiAdj));


% figure;
% subplot(3,1,1)
% plot(x,eta,'r',x,phi,'b',x,phif,'k',x,tauAz,'g',x,tau,'c');%xlim([-100;100])
% subplot(3,1,2)
% plot(x,FluxF0,'r',x,FluxF0+FEv,'--b',x,FluxF,'-.k',x,Fs,'g');%xlim([-100;100])
% subplot(3,1,3)
% plot(x,Lphif0,'r',x,Lphif0+Lev,'b',x,Lphif,'k');%xlim([-100;100])
% %%%%%%%%%%%%%%%%%%%%%%%

    
    shippar.difpot.Amat=0;shippar.difpot.Bvect=0;
    
    if FlagAllSFixed==0
        Nship   = shippar.Nship;
        chiship    =shippar.form.chi;
        chiwaterline=shippar.form.chiWl;
        nu         =shippar.form.nu;%normalized normal vector
        xShip      =shippar.form.xShip;
        Sdraft     =shippar.form.Sdraft;
        Smotion    =shippar.data(:,2);
        
        if Nship==1
            
           
        [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
            shippar.rad.Fpsi,shippar.rad.FpsiEv,~]=...
            funS_radiation_potential_Ndom(shippar,par.x,xShip,Sdraft,0,chiship,chiwaterline,nu,Smotion,-par.bathy',Nship,0,0,0,0);
        
        
%         [psi0,Ma0,dampC0,LpsiEv0,...
%             Fpsi0,FpsiEv0,~]=...
%             funS_radiation_potential(shippar,par.x,xShip,Sdraft,0,chiship,chiwaterline,nu,Smotion,-par.bathy',Nship,0,0,0,0);
%         figure;
%         subplot(3,1,1)
%         plot(par.x,psi0.surge,'b',par.x,shippar.rad.psi.surge,'r')
%         subplot(3,1,2)
%         plot(par.x,psi0.heave,'b',par.x,shippar.rad.psi.heave,'r')
%         subplot(3,1,3)
%         plot(par.x,psi0.pitch,'b',par.x,shippar.rad.psi.pitch,'r')
%         
%         figure;
%         subplot(3,1,1)
%         plot(par.x,Fpsi0.surge,'b',par.x,shippar.rad.Fpsi.surge,'r')
%         subplot(3,1,2)
%         plot(par.x,Fpsi0.heave,'b',par.x,shippar.rad.Fpsi.heave,'r')
%         subplot(3,1,3)
%         plot(par.x,Fpsi0.pitch,'b',par.x,shippar.rad.Fpsi.pitch,'r')
%         
%         
%           [shippar.rad.Ma.x,shippar.rad.Ma.z,shippar.rad.Ma.theta]  
%           [ Ma0.x,Ma0.z,Ma0.theta]
          
          
        elseif Nship==2
            etaWl=[0 0 0 0];dxetaWl=[0 0 0 0];
            [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
                shippar.rad.Fpsi,shippar.rad.FpsiEv,shippar.rad.FpsiS,~]=...
                funS_radiation_potential_2barges(shippar,par.x,xShip,Sdraft,[0 0],chiship,chiwaterline,nu,Smotion,etaWl,dxetaWl,-par.bathy');
        end
    end



% kappa=kappaD./D;
% b=shippar.form.Slength(1)/2;
% indc=closest(x,0);
% indl=closest(x,-b);indr=closest(x,b);
% Ev1=exp(kappa(1)*(x+b));
% Ev1(indc:end)=exp(-kappa(1)*(x(indc:end)-b));
% Ev1(indl+1:indr-1)=0;
% [size(k) size(Ev1)]
% dxEv1=Ifft(1i.*k.*fft(Ev1'));
%
%  figure(100);
%  plot(x,Ev1,'r',x,dxEv1,'b');%,xx,Ev2,'k--',xx,FEv2,'k',xx,u,'c')
%  legend('Ev1','dxEv1');
%  plot_properties;
%

% % inverse operators preparation
% shippar.interpRef=2;
% if shippar.interpRef==2
%      DD=-par.bathy;
%     Om         = str2func(model.dispersion);omAdd      = model.OmFun;
%     Up         = str2func(['Up',model.dispersion(3:end)]);
%
%     [Dmin,Dplus,gammin,gamplus] = SLinv2IP(par.k,Om,omAdd,influx.nu_p,DD,...
%     Proj.Dir,Proj.savename,model.dyn,shippar.form.Sdraft,influx.Hs);
%     [DminCp2,DplusCp2,gamminCp2,gamplusCp2] = SCp2inv2IP(par.k,Om,Up,omAdd,influx.nu_p,DD,...
%     Proj.Dir,Proj.savename,model.dyn,shippar.form.Sdraft,influx.Hs);
%
%     shippar.Oprt.Linvmin         = OpLinv(par.k,Dmin);
%     shippar.Oprt.Linvplus        = OpLinv(par.k,Dplus);
%     shippar.Oprt.Linvgammin      = gammin;
%     shippar.Oprt.Linvgamplus     = gamplus;
%     shippar.Oprt.Cp2invmin       = 1./(Up(par.k,DminCp2,omAdd).^2);
%     shippar.Oprt.Cp2invplus      = 1./(Up(par.k,DplusCp2,omAdd).^2);
%     shippar.Oprt.Cp2invgammin    = gamminCp2;
%     shippar.Oprt.Cp2invgamplus   = gamplusCp2;
%     shippar.Oprt.interp      =2;
%   else
% %    [D_min,D_mid,D_plus,gammin,gammid,gamplus] = Lin3IP(k,Up,Om,omAdd,influx.nu_p,DD,...
% %     Proj.Dir,Proj.savename,model.dyn);
% %     Oprt.Csqmin      = Up(k,D_min,omAdd).^2;
% %     Oprt.Csqmid      = Up(k,D_mid,omAdd).^2;
% %     Oprt.Csqplus     = Up(k,D_plus,omAdd).^2;
% %     Oprt.gammin      = gammin;
% %     Oprt.gammid      = gammid;
% %     Oprt.gamplus     = gamplus;
% %     Oprt.interp      =3;
%  end
%
%
% [shippar.tau_norm]=normalized_tangential_potential(shippar,par,model,bath,Oprt,statusbarObj,jProgressBar);
% shippar.radmethod='minimization';
% shippar.radmethod='interpolation';
% shippar.radmethod='neumannseries';
% if strcmp(shippar.radmethod,'minimization')
%     [shippar.rad]=radiation_potential_minimization(shippar,par,model,bath,Oprt,statusbarObj,jProgressBar);
%     elseif strcmp(shippar.radmethod,'interpolation')
%     shippar.interpRef=2;
%     [shippar.rad]=radiation_potential_interpolation(shippar,par,model,influx,Proj);
% end

% [shippar.rad]=radiation_potential_minimization(shippar,par,model,bath,Oprt,statusbarObj,jProgressBar);


%%% Initial condition
if IVP.type~=1
    shippar.init.phi0=Ifft(IVP.zu0_hat./(1i.*k));
    shippar.init.zeta0       = Ifft(IVP.zeta0_hat);
else
    shippar.init.zeta0       =shippar.form.shapeXcZc(end,:).';%
    shippar.init.phi0        =zeros(size(shippar.init.zeta0));
end


