%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for RunUp  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = RHSShoreShip_etavirt(time,z,model,par,influx,IVP,Oprt,shippar,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig           =influx.gen.INsig;
INsig_skew      =influx.gen.INsig_skew;
gamX            =influx.gen.gamX;
gamX_skew       =influx.gen.gamX_skew;
Indirection     =model.influx.direction;
cfSA            =par.cfSA;
ChiAdj          =influx.gen.nonlinAdj;
jProgressBar    =ProgStatbar.jProgressBar;
statusbarObj    =ProgStatbar.statusbarObj;
aal             =Oprt.aal;
H_minDispInterp =par.interp.H_min1;
H_minShore      =par.interp.H_minShore;
IntCoef         =par.interp.IntCoef;

global flagWarnIVP iterProgBar ITERdt tprev Idstop shipsavevar itersavevar
Visual_ODE='Off';
if IVP.type~=1
    if flagWarnIVP==1
        w = warning('query','last');  %turn off warning odeplot
        warning('off',w.identifier);  %for initial value problem case
        flagWarnIVP=0;
    end
end

if par.IP.check==1
    global dteta dtphihat iterInterior iterInteriordt timeIP uIP etaIP
end

dampchar    = (1-cfSA)';
depth=-par.bathy;
Chi=ChiAdj;

if influx.input.type~=1
    if strcmp(Indirection, 'Uni+')
        Source       = ppval(INsig,time).*gamX - ppval(INsig_skew,time).*gamX_skew;
    elseif strcmp(Indirection, 'Uni-')
        Source       = ppval(INsig,time).*gamX + ppval(INsig_skew,time).*gamX_skew;
    elseif strcmp(Indirection, 'Bi')
        Source       = 2*ppval(INsig,time).*gamX ;
    end
else
    Source=0;
end


%%Progress bar
if time>par.t_init+iterProgBar*par.dt*floor(0.1*par.Nt) %update every 10%
    set(jProgressBar,'Maximum',par.Nt, 'Value',iterProgBar*floor(0.1*par.Nt));
    jProgressBar.setStringPainted( true );
    ETA=remain_time(iterProgBar*floor(0.1*par.Nt),par.Nt);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
    iterProgBar=iterProgBar+1;
end


%%%%%%%%%%%%%%%%%%%
n       = length(k);
z1      = z(1:n);      % etahat
z2      = z(n+1:2*n);  % uhat
eta     = Ifft(z1);
phi       = Ifft(z2);
u_hat   = 1i.*k.*z2;
u       = Ifft(u_hat);

Nship   = shippar.Nship;
sX      = z(2*n+1:2*n+Nship);
Px      = z(2*n+Nship+1:2*n+Nship*2);
sZ      = z(2*n+2*Nship+1:2*n+Nship*3);
Pz      = z(2*n+3*Nship+1:2*n+Nship*4);
sTheta  = z(2*n+4*Nship+1:2*n+Nship*5);
Ptheta  = z(2*n+5*Nship+1:2*n+Nship*6);
Ext     = z(2*n+Nship*6+1:2*n+Nship*6*2);

%%%%%%%%%%%%%%%%%%%%%%%%%%ship param
Nship   = shippar.Nship;
if shippar.calclin==0
[shippar.form]= shapereformsetup(g,x,depth,influx.nu_p,shippar,sX,sZ,sTheta);
end

chiship    =shippar.form.chi;
chiwaterline=shippar.form.chiWl;
dxchiship  =shippar.form.dxchi;
chiWl      =chiwaterline(:,end);
chiS      =chiship(:,end);
chiSWl     =chiS+chiWl;
nutild     =shippar.form.nutild;%unnormalized normal vector
nu         =shippar.form.nu;%normalized normal vector
xShip      =shippar.form.xShip;
Sdraft     =shippar.form.Sdraft;
Slength    =shippar.form.Slength;
chifluid   =(1-chiship(:,end));
chiFS      =chifluid+chiS;
Evmode     =shippar.Evmode;
kappaS     =shippar.form.kappa;


% Amat       =shippar.difpot.Amat;
% Bvect      =shippar.difpot.Bvect;
Smotion    =shippar.data(:,2);

%%%% update operator   
if shippar.calclin==0
    %%% this is applied in calc of ship rad potential
    %%% wave dynamic do not use this
    HtotS=depth+shippar.form.shapeXcZc(end,:).';
    Oprt.interp.p1=ppval(par.interp.IntCoef.Sp_p1,HtotS);
    Oprt.interp.c1=ppval(par.interp.IntCoef.Sp_c1,HtotS);
    Oprt.interp.m1=ppval(par.interp.IntCoef.Sp_m1,HtotS);
    Oprt.interp.p2=ppval(par.interp.IntCoef.Sp_p2,HtotS);
    Oprt.interp.c2=ppval(par.interp.IntCoef.Sp_c2,HtotS);
    Oprt.interp.m2=ppval(par.interp.IntCoef.Sp_m2,HtotS);
    Oprt.interp.p1(HtotS<par.interp.H_min1)=0;
    Oprt.interp.c1(HtotS<par.interp.H_min1)=0;
    Oprt.interp.m1(HtotS<par.interp.H_min1)=0;
    Oprt.interp.p2(HtotS<par.interp.H_min1)=0;
    Oprt.interp.c2(HtotS<par.interp.H_min1)=0;
    Oprt.interp.m2(HtotS<par.interp.H_min1)=0;
    OprtS=Oprt;
    
  %  shippar.tangenflowGlobal=fun_tangential_flow_global_Interp(g,k,Oprt,chiSWl(:,end),'runup');
  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




Upeta   = Ifft(Oprt.Cpeak.*z1);
Upphi     = Ifft(Oprt.Cpeak.*z2);

H=Chi.*eta-par.bathymin;

%%%Define active domain by a characteristic function
WaveChar=ones(size(eta));
WaveChar(H<H_minShore)=zeros;
IndxZerofirst=find(WaveChar==0,1,'first'); 
IndxZerolast=find(WaveChar==0,1,'last');
WaveChar(IndxZerofirst:IndxZerolast)=0;
dampchar(WaveChar==0)=ones;   
if IndxZerofirst==1 %% propagation to the left
dampchar(IndxZerofirst:IndxZerolast+1)=smooth(dampchar(IndxZerofirst:IndxZerolast+1),3);
else
dampchar(IndxZerofirst-1:IndxZerolast)=smooth(dampchar(IndxZerofirst-1:IndxZerolast),3);   
end

% figure;
% plot(x,dampchar);
% pause;

if strcmpi(par.wall.presence,'yes')
    if par.wall.method==1
        Wallchar    =ones(n,1);
    else
        Wallchar    =par.wall.Rho';
    end
    indEndBreak =closest(x,par.wall.position);
else
    Wallchar    =ones(n,1);
    
    if depth(end)==0  %%% overtopping is not covered here.
    indEndBreak =IndxZerofirst;
    else
    indEndBreak =IndxZerolast;
    end
end
ShoreChar=1;%% changed on 12/08/2019 do not use the characteristic function.
HeavChar=Wallchar.*ShoreChar;


%%%% ship potentials
if shippar.calclin==1
    sZInDif=0;
else
    sZInDif=sZ;
end

% [phif,Fs,FEv,Lev,beta,dxi_Kdiff,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3...
%     (shippar.calclin,k,x,xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
%     depth,Nship,Evmode,kappaL,kappaR,Amat,Bvect,nutild,shippar.form.XZrot0(1,:));


if Nship==1
    [phif,tau,Fs,FEv,Lev,beta,dxi_Kdiff,etaWl,dxetaWl,tauAz]=funS_tangential_potential3_Ndom...
        (shippar.calclin,k,x.',shippar.form.xShip,phi,eta,shippar.form.chi,shippar.form.chiWl,...
        -par.bathy,shippar.Evmode,shippar.form.XZrot0,shippar.form.kappa,shippar.tangenflowGlobal,...
        shippar.form.shapeXcZc,shippar.form.shapedat);
    %
    %         [phif,tau,Fs,FEv,Lev,beta,dxi_Kdiff,etaWl,dxetaWl]=funS_tangential_potential3_1...
    %         (shippar.calclin,k,x.',xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
    %         -par.bathy,Evmode,shippar.form.XZrot0,influx.lambda_p,kappaS,...
    %         shippar.tangenflowGlobal,Oprt.L0,Oprt.F0,model.nonlinear,influx.gen.nonlinAdj);
    % %
    % % [phif,Fs,FEv,Lev,beta,dxi_Kdiff,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3...
    %     (shippar.calclin,k,x,xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
    %     -par.bathy',Nship,Evmode,kappaS,Amat,Bvect,nutild,shippar.form.XZrot0(1,:));
    
elseif Nship==2
    [phif,tau,Fs,FEv,Lev,beta,dxi_Kdiff,etaWl,dxetaWl]=funS_tangential_potential3_2barges...
        (shippar.calclin,k,x.',xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
        -par.bathy,Evmode,shippar.form.XZrot0,influx.lambda_p,kappaS,...
        shippar.tangenflowGlobal,Oprt.L0,Oprt.F0,model.nonlinear,influx.gen.nonlinAdj);
end


if shippar.calclin==0
    if any(~strcmp(Smotion(:),'Fixed'))
        %         [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,shippar.rad.Fpsi,shippar.rad.FpsiEv,shippar.rad.FpsiS,dxi_Krad]=...
        %             funS_radiation_potential(shippar,x,xShip,Sdraft,sZ,chiship,chiwaterline,nu,Smotion,depth,Nship,0*etaL,0*etaR,0*dxetaL,0*dxetaR);
        if Nship==1
            [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
                shippar.rad.Fpsi,shippar.rad.FpsiEv,dxi_Krad]=fun_radiation_potential_anzats_Interp(g,k,par.x.',OprtS,shippar.tangenflowGlobal,shippar.form.shapeXcZc(end,:).',shippar.form.dxshapeXcZc(end,:).',Smotion,2,nutild,xShip,shippar.form.XZrot0,chiSWl(:,end),1,0,'runup');
                      
%             [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
%                 shippar.rad.Fpsi,shippar.rad.FpsiEv,dxi_Krad]=...
%                 funS_radiation_potential_Ndom(shippar,par.x,xShip,Sdraft,0,chiship,chiwaterline,nu,Smotion,-par.bathy',Nship,0,0,0,0);
%             
            
            %                 [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
            %                 shippar.rad.Fpsi,shippar.rad.FpsiEv,dxi_Krad]=...
            %                 funS_radiation_potential(shippar,x,xShip,Sdraft,sZ,chiship,...
            %                 chiwaterline,nu,Smotion,-par.bathy',Nship,0,0,0,0);
            
        else
            etaWl=[0 0 0 0];dxetaWl=[0 0 0 0];
            [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
                shippar.rad.Fpsi,shippar.rad.FpsiEv,shippar.rad.FpsiS,dxi_Krad]=...
                funS_radiation_potential_2barges(shippar,par.x,xShip,Sdraft,[0 0],chiship,chiwaterline,nu,Smotion,etaWl,dxetaWl,-par.bathy');
            
        end
        
        
        
    else
        dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
    end
else
    dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
end

[phiRad,Vel,dtetaS,FpsiS,dxi_Krad]=fun_vel_and_radpot_calc(shippar,chiship,chiwaterline,...
         nutild,Pz,Px,Ptheta,beta,Smotion,Nship,dxi_Krad);

phitot=phif+phiRad;


%%%%%%%%%%%%%%%


if strcmp(model.dispersion,'OmSWE')   % Use direct model to avoid interpolation
    if model.nonlinear==2
        Htot        =(depth+Chi.*eta);
        deluH_hat   =fft(Htot.*u.*ShoreChar);
        deletaH_hat = fft(Ifft(g*z1+fft(Chi.*u.^2)./2).*ShoreChar);
        dteta_hat   = -1i*k.*deluH_hat ;
        dtphi_hat   = -deletaH_hat ;
    end    
elseif model.nonlinear>=2
    C2m1 = Oprt.C2m1; C2p1= Oprt.C2p1; C2c1  = Oprt.C2c1;
    C2m2 = Oprt.C2m2; C2p2 = Oprt.C2p2; C2c2 = Oprt.C2c2;
    Cm1 = sqrt(Oprt.C2m1); Cp1= sqrt(Oprt.C2p1); Cc1  = sqrt(Oprt.C2c1);
    Cm2 = sqrt(Oprt.C2m2); Cp2 = sqrt(Oprt.C2p2); Cc2 = sqrt(Oprt.C2c2);
    
    gm1 = IntCoef.m1;   gp1 = IntCoef.p1;  gc1 = IntCoef.c1;
    gm2 = IntCoef.m2;   gp2 = IntCoef.p2;  gc2 = IntCoef.c2;
    
    C_u_hat         = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,u,u_hat);
    C_u_Char        = Ifft(C_u_hat).*HeavChar;C_u_Char_hat=fft(C_u_Char);
    Csq_u_Char_hat  = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,C_u_Char,C_u_Char_hat);
    
    deluH1_hat      = 1/g.*Csq_u_Char_hat;
    gEta=g.*(eta+par.bathyplus);
%     if strcmp(influx.categoryMaxDepth,'Deep water')
%         gEta(H<H_minShore)=g.*H(H<H_minShore);
%     end
    deletaH1_hat    = fft(gEta.*ShoreChar);
    eta(depth<H_minDispInterp)=H(depth<H_minDispInterp);
     
    C2u_hat          = Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,u,u_hat);
    M0u_hat          = -1i.*k.*C2u_hat/g;
    M0u              = Ifft(M0u_hat);
    etaM0uChar       = eta.*M0u.*HeavChar;
    etaM0uChar_hat   = fft(etaM0uChar);
    C2etaM0uChar_hat = Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,etaM0uChar,etaM0uChar_hat);
    M1etaM0uChar_hat = -C2etaM0uChar_hat.*1i.*k./g;
    deluH2_hat       = fft(eta.*u.*HeavChar)+M1etaM0uChar_hat;
    deletaH2_hat     = fft((u.^2 - M0u.^2).*HeavChar)./2;
    
   if model.nonlinear==2
    dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat)));
    dtphi_hat     = -deletaH1_hat -fft(Chi.*Ifft(deletaH2_hat));
    else
        
            M0u_Char             = M0u.*HeavChar;
            dxu                  = Ifft(1i.*k.*u_hat);
            dxetaM0u_hat         = 1i.*k.*fft(eta.*M0u);dxetaMu=Ifft(dxetaM0u_hat);
            CpdxetaM0u_hat       = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,dxetaMu,dxetaM0u_hat);
            CpdxetaM0u_Char      = Ifft(CpdxetaM0u_hat).*HeavChar;
            CpdxetaM0u_Char_hat  = fft(CpdxetaM0u_Char);
            CpCpdxetaM0uChar_hat   = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,CpdxetaM0u_Char,CpdxetaM0u_Char_hat);
            M0dxetaM0uChar_hat     = -1/g.*1i.*k.*CpCpdxetaM0uChar_hat;
            M0dxetaM0uChar         = Ifft(M0dxetaM0uChar_hat);
            etaM0dxetaM0uChar      = eta.*M0dxetaM0uChar;
            etaM0dxetaM0uChar_hat  = fft(etaM0dxetaM0uChar);
            C2etaM0dxetaM0uWallChar_hat=Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,C2m2,gm2,C2p2,gp2,C2c2,gc2,etaM0dxetaM0uChar,etaM0dxetaM0uChar_hat);
            M1etaM0dxetaM0uWallChar_hat=-C2etaM0dxetaM0uWallChar_hat.*1i.*k./g;
            eta2_dxu               = eta.^2.*dxu;   eta2_dxu_hat = fft(eta2_dxu);
            Cp_eta2_dxu_hat        = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,eta2_dxu,eta2_dxu_hat);
            Cp_eta2_dxu_Char       = Ifft(Cp_eta2_dxu_hat).*HeavChar;
            Cp_eta2_dxu_WallChar_hat   = fft(Cp_eta2_dxu_Char);
            CpCpeta2_dxuChar_hat   = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,Cp_eta2_dxu_Char,Cp_eta2_dxu_WallChar_hat);
            M1eta2_dxuChar_hat     = -1/g.*1i.*k.*CpCpeta2_dxuChar_hat;
            deluH3_hat          =  -M1etaM0dxetaM0uWallChar_hat...
                -(1/2).*M1eta2_dxuChar_hat...
                -(1/2).*1i.*k.*fft(eta.^2.*M0u_Char);
            deletaH3_hat        = fft(M0u.*M0dxetaM0uChar)...
                +fft(M0u_Char.*eta.*dxu);
        dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat)));
        dtphi_hat     = -deletaH1_hat -fft(Chi.*Ifft(deletaH2_hat+deletaH3_hat));
    end
end

Upeta(WaveChar==0)=Oprt.Cpeak.*H(WaveChar==0);

Iddamplid=0;
sigma=zeros(size(eta));
if  Nship==2 && Iddamplid==1 
    Xwl=[xShip(1,1) xShip(1,3) xShip(2,1) xShip(2,3)];
    indWl=closest(x,Xwl);
    epsilon=0.05;
    sigma(indWl(2):indWl(3))=epsilon.*influx.nu_p;
end

deluH_hat=deluH1_hat+deluH2_hat;
if model.nonlinear==3
deluH_hat=deluH_hat+deluH3_hat;
end
FluxF=Ifft(deluH_hat);
Lphif_hat=dteta_hat+fft(Lev.*chifluid);
LphiF_hat=(Lphif_hat-fft(sigma.*eta)+fft((FluxF+FEv-Fs-FpsiS).*dxchiship(:,end)));
dteta_hat=  LphiF_hat;

%Breaking process
if strcmp(model.breaking.check,'Yes')
    [B]=breaking_process(eta,u,depth,time,x,par.break,indEndBreak,par.t_init);
    visc=B.*par.break.delb.^2.*(H).*ifft(dteta_hat,'symmetric');
    Flux=-visc.*Ifft(dteta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=zeros(n,1);
    oneperH(H>H_minShore)=1./H(H>H_minShore);
    Rb_hat=fft(oneperH.*dxFlux.*Wallchar)./(1i.*k);
    
    
else
    Rb_hat=0;
end

if par.friction.check==1
    C2m1=par.friction.Cf;%%manning coefficient [0.01;0.04]
    indexf=par.friction.index;
    Rf=zeros(n,1);
    Rf(indexf)=-g*C2m1.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
    Rf(WaveChar==0)=0;
    Rf_hat=fft(Rf)./(1i.*k);
else
    Rf_hat=0;
end

if influx.wind.check==1
    Rw_hat=influx.wind.coef.*z2;
else
    Rw_hat=0;
end

if strcmpi(par.wall.presence,'yes') && par.wall.method==1
    Swall_hat=fun_wallinfluxsource(influx,Oprt,par,eta,time);
else
    Swall_hat=0;
end

dampcharShip=shippar.form.cfSA_ship.';
Sdampcoef=shippar.form.Sdampcoef;
% 
% %%%%%%%%%%%%%%%%%%%%%%% wave dynamic
Cp_phihat   = Op_Cr(Cm1,gm1,Cp1,gp1,Cc1,gc1,Cm2,gm2,Cp2,gp2,Cc2,gc2,phitot,fft(phitot));
etaVirt     = smooth(Ifft(1i.*k.*Cp_phihat)/g,5);

dteta_hat   = (dteta_hat + Swall_hat + Source - fft(Sdampcoef.*dampcharShip.*(eta-etaVirt)+par.dampcoef.*dampchar.*Upeta)).*aal;
dtphi_hat   = (dtphi_hat + Rb_hat +Rf_hat+Rw_hat -fft(Sdampcoef.*dampcharShip.*(phi-phitot)+par.dampcoef.*dampchar.*Upphi)).*aal;

%%%%%%%%%%%%%%%%%%%%%% ship dynamic
if ~any(strcmp(Smotion,'Fixed')) && strcmpi(shippar.moored.check,'yes')
    xi=[sX(1); sZ(1); sTheta(1)];
    if strcmp(Smotion,'Surge')
        Tn=[shippar.moored.Tn1; 0; 0];
    elseif strcmp(Smotion,'Heave')
        Tn=[0; shippar.moored.Tn3; 0];
    elseif strcmp(Smotion,'Pitch')
        Tn=[0; 0; shippar.moored.Tn5];
    else
        Tn=[shippar.moored.Tn1; shippar.moored.Tn3; shippar.moored.Tn5];
    end
    wn2=(2*pi./Tn).^2;
    wn2(Tn==0)=0;
    dtPimoored=((shippar.form.MassS(1).Mat+shippar.rad.Ma.S(1).Mat)*xi).*wn2;
end


if any(strcmp(Smotion,'Surge'))||any(strcmp(Smotion,'Free'))
    dtsX   = Vel.x;
    Zeta   =zeros(Nship,length(x));
    dxZeta =zeros(Nship,length(x));
    dtPx   =zeros(Nship,1);
    
    for ii=1:Nship
        indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
        Zeta(ii,:)   = shippar.form.shapeXcZc(ii,:).*chiship(:,ii).';
        dxZeta(ii,:) = shippar.form.dxshapeXcZc(ii,:);
        dtPx(ii,1)   =0;% -g*trapz(x(indl:indr),Zeta(ii,indl:indr).*nutild.x(ii,indl:indr)).';
        
        if strcmpi(shippar.moored.check,'yes')&& Tn(1)~=0
            dtPx(ii,1)   = dtPx(ii,1)-dtPimoored(1);
        end
       if shippar.form.Cd.X(1)>0
            Cd=shippar.form.Cd.X(1); %Cd=5  ,D=54 and 40, Cd=6  ,D=29; Cd=7  ,D=21;
            dtPx(ii,1)=dtPx(ii,1)-Cd*Sdraft*dtsX*abs(dtsX)/2;
        end
        if shippar.calclin==0
            dtPx(ii,1)   = dtPx(ii,1) -dxi_Kdiff.x(1) -dxi_Krad.x(1);
        end
    end
else
    dtsX   = zeros(Nship,1); dtPx=zeros(Nship,1);%
end

if any(strcmp(Smotion,'Heave'))||any(strcmp(Smotion,'Free'))
    dtsZ   = Vel.z;
    dtPz=zeros(Nship,1);
    
    for ii=1:Nship
        indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
        intchiship  = trapz(x(indl:indr),chiship(indl:indr,ii))';
        dtPz(ii,1)   = -g*(intchiship.*sZ(ii));
        %  dtPz(ii,1)   = -g*Slength(ii).*sZ(ii);
        
        if strcmpi(shippar.moored.check,'yes')&& Tn(2)>0
            dtPz(ii,1)   = -dtPimoored(2);
        end
        
        if shippar.form.Cd.Z(1)>0
            Cd=shippar.form.Cd.Z(1); %Cd=5  ,D=54 and 40, Cd=6  ,D=29; Cd=7  ,D=21;
            dtPz(ii,1)=dtPz(ii,1)-Cd*Sdraft*dtsZ*abs(dtsZ)/2;
        end
        if shippar.calclin==0
            dtPz(ii,1)   =dtPz(ii,1)  -dxi_Kdiff.z(1)-dxi_Krad.z(1);
        end
    end
    
else
    dtsZ   = zeros(Nship,1); dtPz=zeros(Nship,1);%
end

if  any(strcmp(Smotion,'Pitch'))||any(strcmp(Smotion,'Free'))
    dtThe  = Vel.theta;
    Zeta=zeros(Nship,length(x));
    dtPthe=zeros(Nship,1);
    
    M=shippar.form.Mass;
    
    for ii=1:Nship
%         indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
%         Xc=x(closest(x,shippar.form.xShip(ii,2)));
%         Zeta(ii,:)   =shippar.form.shapeXcZc(ii,:).*chiship(:,ii).';
        
        GM=shippar.form.GM; 
        
        dtPthe(ii,1)   = -g*M*GM.*sTheta(ii);
       
        if strcmpi(shippar.moored.check,'yes') && Tn(3)>0
            dtPthe(ii,1)   = -dtPimoored(3) ;
        end
        
        if shippar.form.Cd.Theta(1)>0
            Cd=shippar.form.Cd.Theta(1); 
            dtPthe(ii,1)=dtPthe(ii,1)-Cd*Slength(ii).^4*dtThe*abs(dtThe)/2;
        end
        
        if shippar.calclin==0
            dtPthe(ii,1)   =dtPthe(ii,1) -dxi_Kdiff.theta(1)-dxi_Krad.theta(1);
        end
    end
    
else
    dtThe  = zeros(Nship,1); dtPthe =zeros(Nship,1);
end

dExt  = zeros(Nship*6,1);

dz    = [dteta_hat; dtphi_hat;dtsX;dtPx; dtsZ; dtPz; dtThe; dtPthe; dExt];

if time>=par.t_init+(itersavevar-1)*par.dt
    if any (time<=shipsavevar.time(1:itersavevar-1))
        itersavevar=closest(shipsavevar.time,time);
    end
    shipsavevar.betaz(itersavevar,1:Nship)=beta.z(1:Nship);
    shipsavevar.betax(itersavevar,1:Nship)=beta.x(1:Nship);
    shipsavevar.betatheta(itersavevar,1:Nship)=beta.theta(1:Nship);
    shipsavevar.phiwave(itersavevar,:)=phif;
    
    shipsavevar.dxi_Kdiffx(itersavevar,1:Nship)=dxi_Kdiff.x(1:Nship);
    shipsavevar.dxi_Kdiffz(itersavevar,1:Nship)=dxi_Kdiff.z(1:Nship);
    shipsavevar.dxi_Kdifftheta(itersavevar,1:Nship)=dxi_Kdiff.theta(1:Nship);
    
    if ~strcmpi(shippar.data(1,2),'Fixed')
        shipsavevar.phirad(itersavevar,:)=phiRad;
        for ii=1:Nship
            shipsavevar.addMassS(ii).Mat(itersavevar,:)=reshape(shippar.rad.Ma.S(ii).Mat,1,[]);
            shipsavevar.dampCoefS(ii).Mat(itersavevar,:)=reshape(shippar.rad.dampC.S(ii).Mat,1,[]);
        end
        shipsavevar.dxi_Kradx(itersavevar,1:Nship)=dxi_Krad.x(1:Nship);
        shipsavevar.dxi_Kradz(itersavevar,1:Nship)=dxi_Krad.z(1:Nship);
        shipsavevar.dxi_Kradtheta(itersavevar,1:Nship)=dxi_Krad.theta(1:Nship);
    end
    shipsavevar.time(itersavevar)=time;
    shipsavevar.zVel(itersavevar,1:Nship)=Vel.z;
    shipsavevar.xVel(itersavevar,1:Nship)=Vel.x;
    shipsavevar.thetaVel(itersavevar,1:Nship)=Vel.theta;
    itersavevar=itersavevar+1;
end

if par.IP.check==1
    if time>=par.IP.time(1) && time<=par.IP.time(2)
        if time>=par.IP.time(1)+(iterInteriordt-1)*par.dt*par.IP.time(3)
            dteta(iterInterior,:)=Ifft(dz(1:n));
            dtphihat(iterInterior,:)=dz(n+1:end)./(1i.*k);
            timeIP(iterInterior,1)=time;
            uIP(iterInterior,:)   =u;
            etaIP(iterInterior,:) =eta;
            if any(k==0)
                dtphihat(isnan(dtphihat(iterInterior,:)))=0;
            end
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end

tprev=time;

if strcmp(Visual_ODE,'On')
    if time>ITERdt*par.dt && time>2100
        if ITERdt==1, fig2=figure;end
%         subplot(5,1,1)
%         plot(x,eta,'r',x,par.bathy,'--k',x,ShoreChar,'b:')
%         subplot(5,1,2)
%         plot(x,u,'r')
%         subplot(5,1,3)
%         plot(x,Ifft(deletaH_hat),'r');
%         subplot(5,1,4)
%         plot(x,Ifft(deluH_hat),'r');
%         subplot(5,1,5)
%         plot(x,Ifft(dteta_hat),'r',x,Ifft(dtphi_hat),'b')
          plot(x,shippar.form.shapeXcZc)
        
        title(['t=',num2str(time)]);
       % ylim([-1 1])
        ITERdt=ITERdt+1;
        pause(0.0001);
        
    end
end

