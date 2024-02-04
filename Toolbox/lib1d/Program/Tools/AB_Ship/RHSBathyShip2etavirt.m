%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for Varying-Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dz = RHSBathyShip2etavirt(time,z,model,par,influx,IVP,Oprt,shippar,ProgStatbar)
plot_ode='ye';
%%%%%passing parameter%%%%%%%%%%%%%%
k=par.k; g=par.g; x=par.x;
INsig      =influx.gen.INsig;
INsig_skew =influx.gen.INsig_skew;
gamX       =influx.gen.gamX;
gamX_skew  =influx.gen.gamX_skew;
Indirection=model.influx.direction;
cfSA       =par.cfSA;
ChiAdj     =influx.gen.nonlinAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;

depth    =-par.bathy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar tprev shipsavevar itersavevar Sxprev ITERdt
if IVP.type~=1
    if flagWarnIVP==1
        w = warning('query','last');  %turn off warning odeplot
        warning('off',w.identifier);
        flagWarnIVP=0;
    end
end

if par.IP.check==1
    global dteta dtphihat iterInterior iterInteriordt timeIP uIP etaIP
end


%%Progress bar
if time>par.t_init+iterProgBar*par.dt*floor(0.1*par.Nt) %update every 10%
    set(jProgressBar,'Maximum',par.Nt, 'Value',iterProgBar*floor(0.1*par.Nt));
    jProgressBar.setStringPainted( true );
    ETA=remain_time(iterProgBar*floor(0.1*par.Nt),par.Nt);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
    iterProgBar=iterProgBar+1;
end

dampchar    = (1-cfSA)';

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

%%%%%%%%%%%%%%%%%%%
n       = length(k);
eta_hat = z(1:n);      % etahat
phi_hat = z(n+1:2*n);  % phihat
eta     = Ifft(eta_hat);
phi     = Ifft(phi_hat);
Nship   = shippar.Nship;
sX      = z(2*n+1:2*n+Nship);
Px      = z(2*n+Nship+1:2*n+Nship*2);
sZ      = z(2*n+2*Nship+1:2*n+Nship*3);
Pz      = z(2*n+3*Nship+1:2*n+Nship*4);
sTheta  = z(2*n+4*Nship+1:2*n+Nship*5);
Ptheta  = z(2*n+5*Nship+1:2*n+Nship*6);
Ext     = z(2*n+Nship*6+1:2*n+Nship*6*2);
u_hat   = phi_hat.*1i.*k;
u       = Ifft(u_hat);


if shippar.calclin==0
[shippar.form]= shapereformsetup(g,x,depth,influx.nu_p,shippar,sX,sZ);
end

chiship    =shippar.form.chi;
chiwaterline=shippar.form.chiWl;
dxchiship  =shippar.form.dxchi;
chiWl      =chiwaterline(:,end);
chiSl      =chiship(:,end);
chiSWl     =chiSl+chiWl;
nutild     =shippar.form.nutild;%unnormalized normal vector
nu         =shippar.form.nu;%normalized normal vector
xShip      =shippar.form.xShip;
Sdraft     =shippar.form.Sdraft;
Slength    =shippar.form.Slength;
chifluid   =(1-chiship(:,end));
chiFS      =chifluid+chiSl;
Evmode     =shippar.Evmode;
kappaS     =shippar.form.kappa;

% Amat       =shippar.difpot.Amat;
% Bvect      =shippar.difpot.Bvect;
Smotion    =shippar.data(:,2);

%%%% update operator   
if shippar.calclin==0
    HtotS=depth+shippar.form.shapeXcZc(end,:).';
    Oprt.S.gammin=ppval(Oprt.S.Sp_gammin,HtotS);
    Oprt.S.gamplus=ppval(Oprt.S.Sp_gamplus,HtotS);
    Oprt.S.gammid=ppval(Oprt.S.Sp_gammid,HtotS);
  %  shippar.tangenflowGlobal=fun_tangential_flow_global_Interp(g,k,Oprt,par.cfSA.','Flat');
end

%%%%%%%%%%%%%%%%%%%%%%% ship part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                shippar.rad.Fpsi,shippar.rad.FpsiEv,dxi_Krad]=fun_radiation_potential_anzats_Interp(g,k,par.x.',Oprt.S,shippar.tangenflowGlobal,shippar.form.shapeXcZc(end,:).',shippar.form.dxshapeXcZc(end,:).',Smotion,2,nutild,xShip,shippar.form.XZrot0,chiSWl(:,end),1,0,'B');

                
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

%%%%%%%%%%%%%%%%%%%%%%% wave part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model.evol,'HS')
    [Lphif_hat,FluxPhi_hat]=funOprt_Lbathy(g,k,Oprt,eta,u,model.nonlinear,ChiAdj);
    
    deletaH_hat    = g.*eta_hat;
    if model.nonlinear~=1
        M0u_hat      = Op_M0(k,g,Oprt,u_hat);
        M0u          = Ifft(M0u_hat);
        deletaH_hat  =deletaH_hat+ fft((u.^2 - M0u.^2))./2;
        if model.nonlinear~=2
        disp('third order is not available yet.') 
        end
    end
end

Iddamplid=0;
sigma=zeros(size(eta));
if  Nship==2 && Iddamplid==1 
    Xwl=[xShip(1,1) xShip(1,3) xShip(2,1) xShip(2,3)];
    indWl=closest(x,Xwl);
    epsilon=0.05;
    sigma(indWl(2):indWl(3))=epsilon.*influx.nu_p;
end
%%%%%%%%%%%%%%%%%%%%%%% ship part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FluxF=Ifft( FluxPhi_hat);

Lphif_hat=Lphif_hat+fft(Lev.*chifluid);
LphiF_hat=(Lphif_hat-fft(sigma.*eta)+fft((FluxF+FEv-Fs-FpsiS).*dxchiship(:,end)));

dteta_hat       =  LphiF_hat;
dtphi_hat       =  -deletaH_hat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if strcmp(model.breaking.check,'Yes')
    
     H=(depth+eta);
    if strcmpi(par.wall.presence,'yes')
        indEndBreak =closest(x,par.wall.position);
    else
        indEndBreak=[];
    end
    [B] = breaking_process(eta,u,depth,time,x,par.break,indEndBreak,par.t_init);
    
    visc=B.*par.break.delb.^2.*(H).*ifft(dteta_hat,'symmetric');
    Flux=-visc.*Ifft(dteta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=1./H;
    Rb_hat=fft(oneperH.*dxFlux)./(1i.*k);
else
    Rb_hat=0;
end

if par.friction.check==1
     H=depth+eta;
    Cm=par.friction.Cf;%%manning coefficient [0.01;0.04]
    indexf=par.friction.index;
    Rf=zeros(n,1);
    Rf(indexf)=-g*Cm.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
    Rf_hat=fft(Rf)./(1i.*k);
    Rf_hat(isnan(Rf_hat))=0;
else
    Rf_hat=0;
end

if influx.wind.check==1
    if time>=influx.wind.tinterv(1)&&time<=influx.wind.tinterv(2)
        Rw=zeros(n,1);
        Rw_temp=influx.wind.coef.*Ifft(Oprt.Omd.*z2);
        indexXi=influx.wind.indexX(1);indexXf=influx.wind.indexX(2);
        Rw(indexXi:indexXf)=Rw_temp(indexXi:indexXf);
        Rw_hat=fft(Rw);
    else
        Rw_hat=0;
    end
else
    Rw_hat=0;
end

if strcmpi(par.wall.presence,'yes')
    Swall_hat=fun_wallinfluxsource(influx,par,eta,time);
else
    Swall_hat=0;
end

dampcharShip=shippar.form.cfSA_ship.';
Sdampcoef=shippar.form.Sdampcoef;
etaVirt=smooth(Ifft(1i.*OpOminterp(Oprt,fft(phitot)))/g,5);
dteta_hat   = (dteta_hat + Swall_hat + Source - fft((Sdampcoef.*dampcharShip.*(eta-etaVirt)+par.dampcoef.*dampchar.*eta))).*aal;
dtphi_hat   = (dtphi_hat + Rb_hat +Rf_hat+Rw_hat -  fft((Sdampcoef.*dampcharShip.*(phi-phitot)+par.dampcoef.*dampchar.*phi))).*aal;


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
        if strcmpi(shippar.moored.check,'yes')&& Tn(2)~=0
            dtPz(ii,1)   = -dtPimoored(2);
        else
            indl=closest(x,shippar.form.xShip(ii,1));indr=closest(x,shippar.form.xShip(ii,3));
            intchiship  = trapz(x(indl:indr),chiship(indl:indr,ii))';
            dtPz(ii,1)   = -g*(intchiship.*sZ(ii));
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
       
        %dtPthe(ii,1)   = -g*trapz(x(indl:indr)-Xc,Zeta(ii,indl:indr).*nutild.theta(ii,indl:indr)).';%-g.*trapz(x(indl:indr)-Xc,Zeta(ii,indl:indr).*(x(indl:indr)-Xc));%
        %         if shippar.calclin==0
        %         dtPthe(ii,1)   =dtPthe(ii,1) -dxi_Kdiff.theta(1)-dxi_Krad.theta(1);
        %         end
        if shippar.form.Cd.Theta(1)>0
            Cd=shippar.form.Cd.Theta(1); %Cd=5  ,D=54 and 40, Cd=6  ,D=29; Cd=7  ,D=21;
            dtPthe(ii,1)=dtPthe(ii,1)-Cd*Sdraft*dtThe*abs(dtThe)/2;
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

if strcmp(plot_ode,'yes')
    if time>ITERdt*par.dt  && time>20
        if ITERdt==1, fig2=figure;end
        plot(x,phi,'b',x,phitot,'--r',x,dphidtS,'g');
        pause;
        ITERdt=ITERdt+1;
    end
end

