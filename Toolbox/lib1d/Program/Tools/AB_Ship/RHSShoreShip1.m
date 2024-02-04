%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for RunUp  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = RHSShoreShip1(time,z,model,par,influx,IVP,Oprt,shippar,ProgStatbar)
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig       =influx.gen.INsig;
INsig_skew  =influx.gen.INsig_skew;
gamX        =influx.gen.gamX;
gamX_skew   =influx.gen.gamX_skew;
Indirection =model.influx.direction;
cfSA        =par.cfSA;
ChiAdj      =influx.gen.nonlinAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;
H_min      =par.interp.H_min1;
H_minShore  =par.interp.H_minShore;
IntCoef     =par.interp.IntCoef;

global flagWarnIVP iterProgBar ITERdt tprev itersavevar shipsavevar
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

dampchar    = (1-cfSA).';
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
u       = Ifft(z2);

phihat =z2./(1i.*k);
phi    = Ifft(phihat);
phiprev=phi;

Nship   = shippar.Nship;
sX      = z(2*n+1:2*n+Nship);
Px      = z(2*n+Nship+1:2*n+Nship*2);
sZ      = z(2*n+2*Nship+1:2*n+Nship*3);
Pz      = z(2*n+3*Nship+1:2*n+Nship*4);
sTheta  = z(2*n+4*Nship+1:2*n+Nship*5);
Ptheta  = z(2*n+5*Nship+1:2*n+Nship*6);
Ext     = z(2*n+Nship*6+1:2*n+Nship*6*2);

chiship    =shippar.form.chi;
chiwaterline=shippar.form.chiWl;
dxchiship  =shippar.form.dxchi;
chiWl      =chiwaterline(:,end);
chiSl      =chiship(:,end);
nutild     =shippar.form.nutild;%unnormalized normal vector
nu         =shippar.form.nu;%normalized normal vector
xShip      =shippar.form.xShip;
Sdraft     =shippar.form.Sdraft;
Slength    =shippar.form.Slength;
chifluid   =(1-chiship(:,end));
chiFS      =chifluid+chiSl;
Evmode     =shippar.Evmode;
kappaL      =shippar.difpot.kappaL; 
kappaR      =shippar.difpot.kappaR;

Amat       =shippar.difpot.Amat;
Bvect      =shippar.difpot.Bvect;
Smotion    =shippar.data(:,2);


%%%%%%%%%%%%%%%%%%%%%%% ship part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if shippar.calclin==1
    sZInDif=0;
else
    sZInDif=sZ;
end

[phif,phitan,phitanZ0,Fs,FEv,Lev,beta,dxi_Kdiff,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3...
    (shippar.calclin,Oprt,k,x,xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
    depth.',Nship,Evmode,kappaL,kappaR,time,Amat,Bvect,Smotion,nutild,shippar.form.XZrot0(1,:));

phif=phif.*chifluid+phitan.*(1-chifluid);

if shippar.calclin==0
    if any(~strcmp(Smotion(:),'Fixed'))
        [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,shippar.rad.Fpsi,shippar.rad.FpsiEv,shippar.rad.FpsiS,dxi_Krad]=...
            funS_radiation_potential(shippar,x,xShip,Sdraft,sZ,chiship,chiwaterline,nu,Smotion,depth.',Nship,etaL,etaR,dxetaL,dxetaR);
    else
        dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
        dxi_Kdiff.z(1)=0; dxi_Kdiff.x(1)=0;dxi_Kdiff.theta(1)=0;
    end
else
    dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
end

[phiRad,Vel,dtetaS,LpsiEv,FpsiEv,FpsiS,PsiDtNPsi,dxi_Krad]=fun_vel_and_radpot_calc(k,Oprt.L0,eta,x,shippar,phif,phitan,chiship,...
    chiWl,dxchiship,nu,nutild,Pz,Px,Ptheta,model.nonlinear,depth,Sdraft,sZ,Slength,beta,Smotion,Nship,dxi_Krad);

phitot=phif+phiRad;
%dphidtS=(phitot-phiprev)/par.dt;
utot=Ifft(1i.*k.*fft(phitot));
dudtS=(utot-u)/par.dt;


H=Chi.*eta+depth;
Upeta   = Ifft(Oprt.Cpeak.*z1);
Upu     = Ifft(Oprt.Cpeak.*z2);


HeavH=Heaviside(H-H_minShore);
indDI=find(HeavH==0,1,'first');
indDF=find(HeavH==0,1,'last');
ShoreChar=HeavH;

if isempty(indDI)
    indWI=1; indWF=n;
    indEndBreak=n;
    FlagOvertop=0;
    signProp=0;
elseif indDI==1  % propagation direction to the left
    ShoreChar(indDI:indDF)=zeros;
    indWI=indDF+1; indWF=n;
    indEndBreak=indWI-10;
    FlagOvertop=0;
    signProp=-1;
elseif indDF==n
%     HeavHC2=Heaviside(H(indDI+1:indDF)-H_minShore);
%     indDI=indDI+find(HeavHC2==0,1,'first');
    
    ShoreChar(indDI:indDF)=zeros;
    indWI=1; indWF=indDI-1;
    indEndBreak=indWF+10;
    FlagOvertop=0;
    signProp=1;
elseif indDI~=1 && indDF~=n
    %HeavH(H<50*H_min)=0;
    indW=find(HeavH==1);
    indD=find(HeavH<H_minShore);
    ShoreChar=HeavH;
    indEndBreak=n;
    FlagOvertop=1;
    signProp=0;
end

if isfield(par,'wall')
    if strcmp(par.wall.presence,'Yes')
        Wallchar    =par.wall.Rho';
        indEndBreak =closest(x,par.wall.position);
    else
        Wallchar    =ones(n,1);
    end
end

ShoreChar=1;%% changed on 12/08/2019 do not use the characteristic function.
indWI=1;indWF=n;
HeavChar=Wallchar.*ShoreChar;

if strcmp(model.dispersion,'OmSWE')   % Use direct model to avoid interpolation
    if model.nonlinear==2
        Htot        =(depth+Chi.*eta);
        deluH_hat   =fft(Htot.*Ifft(z2).*ShoreChar);
        deletaH_hat = fft(Ifft(g*z1+fft(Chi.*u.^2)./2).*ShoreChar);
        dteta_hat     = -1i*k.*deluH_hat ;
        dtu_hat       = -1i*k.*deletaH_hat ;
    end
elseif model.nonlinear==2
    C2m1         = Oprt.C2m1; 
    C2p1         = Oprt.C2p1; 
    C2c1         = Oprt.C2c1; 
    C2m2         = Oprt.C2m2; 
    C2p2         = Oprt.C2p2; 
    C2c2         = Oprt.C2c2; 
    
    gm1=zeros(n,1);gp1=gm1;gc1=gp1;
    gm2=zeros(n,1);gp2=gm2;gc2=gp2;
    
    
    if FlagOvertop==0
        gm1(indWI:indWF)=ppval(IntCoef.m1,H(indWI:indWF));
        gp1(indWI:indWF)=ppval(IntCoef.p1,H(indWI:indWF));
        gc1(indWI:indWF)=ppval(IntCoef.c1,H(indWI:indWF));
        gm2(indWI:indWF)=ppval(IntCoef.m2,H(indWI:indWF));
        gp2(indWI:indWF)=ppval(IntCoef.p2,H(indWI:indWF));
        gc2(indWI:indWF)=ppval(IntCoef.c2,H(indWI:indWF));
    else
        gm1(indW)=ppval(IntCoef.m,H(indW));
        gp1(indW)=ppval(IntCoef.p,H(indW));
        gc1(indW)=ppval(IntCoef.c,H(indW));
    end
    
   
    CCu_hat= Op_C2r(C2m1,gm1,C2p1,gp1,C2c1,gc1,...
                   C2m2,gm2,C2p2,gp2,C2c2,gc2,u,z2);
    CCu_ShoreChar    =Ifft(CCu_hat).*ShoreChar.*Wallchar;
    uChar           =u.*ShoreChar.*Wallchar;
    CCu_ShoreChar(H<H_min) = g*H(H<H_min).*uChar(H<H_min);
    CCu_ShoreChar_hat=fft(CCu_ShoreChar);
 
    kC2_u=Ifft(k.*CCu_ShoreChar_hat);
    
    deluH_hat   =  1/g.*CCu_ShoreChar_hat;
    deletaH_hat =  g.*z1+ fft(Ifft(fft(u.^2/2-(kC2_u./g).^2/2)).*Chi.*Wallchar);
    deletaH_hat = fft(Ifft(deletaH_hat).*ShoreChar);
    
    
      %%%%%%%%%%%%%%%%%%%%%%% wave-ship part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   FluxF=Ifft( deluH_hat);
   LphiF_hat=fft(Ifft((-1i.*k.*(deluH_hat)+fft(Lev)+fft((FluxF+FEv-Fs-FpsiS).*dxchiship(:,end))).*aal).*chifluid+dtetaS);
    
    dteta_hat = LphiF_hat;     
    dtu_hat   = -fft(Ifft(1i.*k.*deletaH_hat).*chifluid-0*dudtS.*(1-chifluid)).*aal ;
        
elseif model.nonlinear>2
    Cm2= Oprt.Cm.^2;  Cp2=Oprt.Cp.^2;  Cc2=Oprt.Cc.^2;
    C2m1= Oprt.Cm;  C2p1=Oprt.Cp;  C2c1=Oprt.Cc;
    gm1 = IntCoef.m;   gp1 = IntCoef.p;  gc1 = IntCoef.c;
    C_u_hat         = Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,u,z2);
    C_u_Char        = Ifft(C_u_hat).*HeavChar;C_u_Char_hat=fft(C_u_Char);
    Csq_u_Char_hat  = Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,C_u_Char,C_u_Char_hat);
    
    
    deluH1_hat      = 1/g.*Csq_u_Char_hat;
    deletaH1_hat    = g.*fft(eta.*ShoreChar);
    
    
    M0u_hat          = Op_M03p(k,g,Cm2,gm1,Cp2,gp1,Cc2,gc1,z2);
    M0u              = Ifft(M0u_hat);
    etaM0uChar_hat   = fft(eta.*M0u.*HeavChar);
    M1etaM0uChar_hat = Op_M13p(k,g,Cm2,gm1,Cp2,gp1,Cc2,gc1,etaM0uChar_hat);
    deluH2_hat       = fft(eta.*u.*HeavChar)+M1etaM0uChar_hat;
    deletaH2_hat     = fft((u.^2 - M0u.^2).*HeavChar)./2;
    
    if model.nonlinear==2
        dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*deluH2_hat));
        dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*deletaH2_hat));
    else
        M0u_Char             = M0u.*HeavChar;
        dxu                  = Ifft(1i.*k.*z2);
        dxetaM0u_hat         = 1i.*k.*fft(eta.*M0u);dxetaMu=Ifft(dxetaM0u_hat);
        CpdxetaM0u_hat       = Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,dxetaMu,dxetaM0u_hat);
        CpdxetaM0u_Char      = Ifft(CpdxetaM0u_hat).*HeavChar;
        CpdxetaM0u_Char_hat  = fft(CpdxetaM0u_Char);
        M0dxetaM0uChar_hat     = -1/g.*1i.*k.*Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,CpdxetaM0u_Char,CpdxetaM0u_Char_hat);
        M0dxetaM0uChar         = Ifft(M0dxetaM0uChar_hat);
        etaM0dxetaM0uChar_hat  = fft(eta.*M0dxetaM0uChar);
        M1etaM0dxetaM0uWallChar_hat= Op_M13p(k,g,Cm2,gm1,Cp2,gp1,Cc2,gc1,etaM0dxetaM0uChar_hat);
        eta2_dxu               = eta.^2.*dxu;   eta2_dxu_hat = fft(eta2_dxu);
        Cp_eta2_dxu_hat        = Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,eta2_dxu,eta2_dxu_hat);
        Cp_eta2_dxu_Char       = Ifft(Cp_eta2_dxu_hat).*HeavChar;
        Cp_eta2_dxu_WallChar_hat   = fft(Cp_eta2_dxu_Char);
        M1eta2_dxuChar_hat     = -1/g.*1i.*k.*Op_C(C2m1,gm1,C2p1,gp1,C2c1,gc1,Cp_eta2_dxu_Char,Cp_eta2_dxu_WallChar_hat);
        deluH3_hat          =  M1etaM0dxetaM0uWallChar_hat...
            -(1/2).*M1eta2_dxuChar_hat...
            -(1/2).*1i.*k.*fft(eta.^2.*M0u_Char);
        deletaH3_hat        = fft(M0u.*M0dxetaM0uChar)...
            +fft(M0u_Char.*eta.*dxu);
        
        if model.nonlinear==3
            dteta_hat   = -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat)));
            dtu_hat     = -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat)));
        else
            %             dx2eta_hat                  =1i.*k.*1i.*k.*z1;
            %             dx2eta                      =Ifft(dx2eta_hat);
            %             dxetaM0u_hat                =1i.*k.*fft(eta.*M0u);
            %             M0dxetaM0u_hat              =Op_M0(k,g,Cm2,gm,Cp2,gp,dxetaM0u_hat);
            %             M0dxetaM0u                  =Ifft(M0dxetaM0u_hat);
            %             dxetaM0dxetaM0uChar_hat     =1i.*k.*fft(eta.*M0dxetaM0u.*HeavChar);
            %             M0dxetaM0dxetaM0uChar_hat   =Op_M0(k,g,Cm2,gm,Cp2,gp,dxetaM0dxetaM0uChar_hat);
            %             M0dxetaM0dxetaM0uChar       =Ifft(M0dxetaM0dxetaM0uChar_hat);
            %             etaM0dxetaM0dxetaM0uChar_hat=fft(eta.*M0dxetaM0dxetaM0uChar);
            %             M1etaM0dxetaM0dxetaM0uChar_hat=Op_M1(k,g,Cm2,gm,Cp2,gp,etaM0dxetaM0dxetaM0uChar_hat);
            %
            %             dxeta2dxuChar_hat         =1i.*k.*fft(eta.^2.*dxu.*HeavChar);
            %             M0dxeta2dxuChar_hat       =Op_M0(k,g,Cm2,gm,Cp2,gp,dxeta2dxuChar_hat);
            %             M0dxeta2dxuChar           =Ifft(M0dxeta2dxuChar_hat);
            %             etaM0dxeta2dxuChar_hat    =fft(eta.*M0dxeta2dxuChar);
            %             M1etaM0dxeta2dxuChar_hat  =Op_M1(k,g,Cm2,gm,Cp2,gp,etaM0dxeta2dxuChar_hat);
            %
            %             dxeta3dxM0uChar_hat      =1i.*k.*fft(eta.^3.*HeavChar.*Ifft(1i.*k.*M0u_hat));
            %             M1dxeta3dxM0uChar_hat    =Op_M1(k,g,Cm2,gm,Cp2,gp,dxeta3dxM0uChar_hat);
            %
            %             eta2dx2etaM0uChar_hat    =fft(HeavChar.*eta.^2.*dx2eta.*M0u);
            %             M1eta2dx2etaM0uWallChar_hat  =Op_M(k,g,Cm2,gm,Cp2,gp,eta2dx2etaM0uChar_hat);
            %
            %             deluH4_hat     = (1/3).*1i.*k.*fft(eta.^3.*dxu.*HeavChar)...
            %                 +M1etaM0dxetaM0dxetaM0uChar_hat...
            %                 +(1/2).*1i.*k.*fft(HeavChar.*eta.^2.*M0dxetaM0u)...
            %                 +(1/2).*M1etaM0dxeta2dxuChar_hat...
            %                 +(1/3).*M1dxeta3dxM0uChar_hat...
            %                 +(1/2).*M1eta2dx2etaM0uWallChar_hat;
            %             deletaH4_hat     = (1/2).*fft( -(eta.^2.*dxu.^2.*HeavChar)...
            %                 -HeavChar.*M0dxetaM0u.^2 ...
            %                 -2.*(M0u.*M0dxetaM0dxetaM0uChar)...
            %                 -2.*(eta.*dxu.*M0dxetaM0u.*HeavChar)...
            %                 -(M0u.*M0dxeta2dxuChar)...
            %                 +(eta.^2.*Ifft(1i.*k.*M0u_hat).^2.*HeavChar)...
            %                 -(eta.*dx2eta.*M0u.^2.*HeavChar))...
            %                 -(1/4).*1i.*k.*1i.*k.*fft(eta.^2.*M0u.^2.*HeavChar);
            %
            %             dteta_hat=  -1i*k.*deluH1_hat -fft(Chi.*Ifft(1i.*k.*(deluH2_hat+deluH3_hat+deluH4_hat)));
            %             dtu_hat  =   -1i*k.*deletaH1_hat -fft(Chi.*Ifft(1i.*k.*(deletaH2_hat+deletaH3_hat+deletaH4_hat)));
        end
    end
end
if ~isempty(indDI)
    if FlagOvertop==0
        dampchar(indDI:indDF)=ones;
        Upeta(indDI:indDF)   =H(indDI:indDF);
        Upu (indDI:indDF)    =u(indDI:indDF);
    else
        dampchar(indD)=ones;
        Upeta(indD)   =H(indD);
        Upu (indD)    =u(indD);
    end
end


%Breaking process
if strcmp(model.breaking.check,'Yes')
    [B]=breaking_process(eta,u,depth,time,x,par.break,indEndBreak);
    visc=B.*par.break.delb.^2.*(H).*ifft(dteta_hat,'symmetric');
    Flux=-visc.*Ifft(dteta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=zeros(n,1);
    oneperH(H>H_minShore)=1./H(H>H_minShore);
    Rb_hat=fft(oneperH.*dxFlux.*Wallchar);
    %         if time>30
    %         subplot(4,1,2)
    %         Rb=Ifft(Rbhat);
    %         plot(x,Rb,'c',x,visc,'--r',x,oneperH,':b',x,u,'k',x,ShoreChar,'g--')
    %         title(['time= ', num2str(time)])
    %         xlim([20;100]);ylim([-10;20]);
    %         pause(0.001)
    %         end
    
    
else
    Rb_hat=0;
end

if par.friction.check==1
    C2m1=par.friction.Cf;%%manning coefficient [0.01;0.04]
    indexf=par.friction.index;
    Rf=zeros(n,1);
    Rf(indexf)=-g*C2m1.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
    Rf(HeavH==0)=0;
    Rf_hat=fft(Rf);
else
    Rf_hat=0;
end

if influx.wind.check==1
    Rw_hat=influx.wind.coef.*z2;
else
    Rw_hat=0;
end


dteta_hat = dteta_hat +( Source -(1+signProp*1i.*k).*fft(Upeta.*dampchar)).*aal;
dtu_hat   = dtu_hat+(Rb_hat+Rf_hat+Rw_hat-(1+signProp*1i.*k).*fft(Upu.*dampchar)).*aal;


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
    dtPimoored=((shippar.form.MassMat)*xi).*(2*pi./Tn).^2;
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
        dtPx(ii,1)   = -g*trapz(x(indl:indr),Zeta(ii,indl:indr).*nutild.x(ii,indl:indr)).'; 
        
        if strcmpi(shippar.moored.check,'yes')
            dtPx(ii,1)   = dtPx(ii,1)-dtPimoored(1);
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
        if strcmpi(shippar.moored.check,'yes')
            dtPz(ii,1)   = -dtPimoored(2);
        else
            indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
            intchiship  = trapz(x(indl:indr),chiship(indl:indr,ii))';
            dtPz(ii,1)   = -g*(intchiship.*sZ(ii));
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
    for ii=1:Nship
        indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
        Xc=x(closest(x,shippar.form.xShip(ii,2)));
        Zeta(ii,:)   =shippar.form.shapeXcZc(ii,:).*chiship(:,ii).';
        
        dtPthe(ii,1)   = -g*trapz(x(indl:indr)-Xc,Zeta(ii,indl:indr).*nutild.theta(ii,indl:indr)).';%-g.*trapz(x(indl:indr)-Xc,Zeta(ii,indl:indr).*(x(indl:indr)-Xc));%
                  if shippar.calclin==0
                     dtPthe(ii,1)   =dtPthe(ii,1) -dxi_Kdiff.theta(1)-dxi_Krad.theta(1);
                  end
        if strcmpi(shippar.moored.check,'yes')
            dtPthe(ii,1)   = dtPthe(ii,1)-dtPimoored(3) ;
        end
        
        
    end
    
else
    dtThe  = zeros(Nship,1); dtPthe =zeros(Nship,1);
end

dExt  = zeros(Nship*6,1);

dz    = [dteta_hat; dtu_hat;dtsX;dtPx; dtsZ; dtPz; dtThe; dtPthe; dExt];

if time>=par.t_init+(itersavevar-1)*par.dt
    if any (time<=shipsavevar.time(1:itersavevar-1))
        itersavevar=closest(shipsavevar.time,time);
    end
    shipsavevar.betaz(itersavevar,1:Nship)=beta.z(1:Nship);
    shipsavevar.betax(itersavevar,1:Nship)=beta.x(1:Nship);
    shipsavevar.betatheta(itersavevar,1:Nship)=beta.theta(1:Nship);
    shipsavevar.phiwave(itersavevar,:)=phif;
    if ~strcmpi(shippar.data(1,2),'Fixed')
        shipsavevar.phirad(itersavevar,:)=phiRad;
        shipsavevar.addMass(itersavevar,:)=reshape(shippar.rad.Ma.Mat,1,[]);
        shipsavevar.dampCoef(itersavevar,:)=reshape(shippar.rad.dampC,1,[]);
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
    if time>ITERdt*par.dt
        %if ITERdt==1, fig2=figure;end;
        plot(x,eta,'r',x,par.bathy,'--k',x,ShoreChar,'b:')
        title(['t=',num2str(time)]);
        % ylim([-0.2 0.2])
        ITERdt=ITERdt+1;
        pause(0.0001);
        
    end
end

