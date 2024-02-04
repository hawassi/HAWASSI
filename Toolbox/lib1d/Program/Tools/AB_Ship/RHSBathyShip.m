%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for Varying-Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function dz = RHSBathyShip(time,z,model,par,influx,IVP,Oprt,shippar,ProgStatbar)
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
depth      =par.depth;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
aal         =Oprt.aal;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar ITERdt shipsavevar itersavevar
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
end;


Nship   = shippar.Nship;
n       = length(k);
eta_hat = z(1:n);      % etahat
phi_hat = z(n+1:2*n);  % uhat
eta     = Ifft(eta_hat);
phi     = Ifft(phi_hat);
sX      = z(2*n+1:2*n+Nship);
Px      = z(2*n+Nship+1:2*n+Nship*2);
sZ      = z(2*n+2*Nship+1:2*n+Nship*3);
Pz      = z(2*n+3*Nship+1:2*n+Nship*4);
The     = z(2*n+4*Nship+1:2*n+Nship*5);
PThe    = z(2*n+5*Nship+1:2*n+Nship*6);
Ext     = z(2*n+Nship*6+1:2*n+Nship*6*2);

 phiprev=phi;
   %etaprev=eta;
%  indlprev=closest(x,Sxprev(1,1));indrprev=closest(x,Sxprev(1,3));
%  xS=sX(1)+shippar.form.xShip(1,2);
%  [shippar.form]= shapeformsetup(x,shippar,xS,sZ(1));

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
kappa      =shippar.difpot.kappa;
Amat       =shippar.difpot.Amat;
Bvect      =shippar.difpot.Bvect;
Smotion    =shippar.data(:,2);

if strcmp(model.evol,'HS')
   
    %%%%% phi fluid whole domain
    DtNphi=Ifft(funOprtBathy_L(g,k,Oprt,eta,phi,model.nonlinear,ChiAdj)).*(1-chiship(:,end));
    phif_hat=funOprtBathy_Linv(g,k,Oprt,eta,DtNphi,model.nonlinear).*aal;%Oprt.aalloc;
    phif=Ifft(phif_hat);
    
    [phiRad,beta,Vel,dtetaS]=fun_vel_momenta_and_radpot_calc(x,shippar,phif,chiship,nu,Pz,Px);
   
    phitot=phif+phiRad;
    dphidtS=((phitot-phi)/par.dt);
  

    %delphiH_hat  = fft(DtNphi+dtetaS);
    delphiH_hat  =funOprtBathy_L(g,k,Oprt,eta,phitot,model.nonlinear,ChiAdj);
    deletaH_hat  = fft(g*eta.*(1-chiship(:,end))-dphidtS.*chiship(:,end));
    if model.nonlinear~=1
       phitot_hat=fft(phitot);
       dxphi=Ifft(1i.*k.*phitot_hat);
       L0phi=Ifft(funOprtBathy_L0(g,k,Oprt,phitot)); 
       deletaH_hat  = deletaH_hat+0.5.*fft((dxphi.^2-L0phi.^2).*(1-chiship(:,end)).*ChiAdj);  
       if model.nonlinear~=2
       L0etaL0phi_hat=funOprtBathy_L0(g,k,Oprt,eta.*L0phi);    
       deletaH_hat  = deletaH_hat+ fft((L0phi.*(eta.*Ifft(1i.*k.*fft(dxphi))+...
                     Ifft(L0etaL0phi_hat))).*(1-chiship(:,end)).*ChiAdj);   
       end
    end
   
    dteta_hat     = delphiH_hat ;
    dtphi_hat     = -deletaH_hat ;
end


dteta_hat   = (dteta_hat + Source - fft(par.dampcoef.*eta.*dampchar)).*aal;
dtphi_hat     = (dtphi_hat  - fft(par.dampcoef.*phi.*dampchar)).*aal;

if any(strcmp(shippar.data(:,2),'Heave'))
dtsZ   = Vel.z;
% intchiship  = trapz(x,chiship(:,1:Nship))';
% dtPz   = -g*(intchiship.*sZ);

ZetaS=shippar.reform.shapeXcZc';
dxZetaS=shippar.reform.dxshapeXcZc';
Mass=shippar.form.Mass;
nuz    =nu.z(1:Nship,:)';
nuz_bar=nuz./sqrt(1+dxZetaS(:,1:Nship).^2);   
dtPz   =-g*Mass-g.*trapz(x,ZetaS(:,1:Nship).*nuz_bar);
   
    if model.nonlinear~=1
            Lphitot=Ifft(funOprtBathy_L(g,k,Oprt,eta,phitot,model.nonlinear,ChiAdj));
            Ks=0.5*Lphitot.*phitot;
        dtPz   =dtPz-trapz(x,repmat(Ks,1,Nship).*nu.z(1:Nship,:)')';
    end    
else
    dtsZ   = zeros(Nship,1); dtPz   = zeros(Nship,1);
end
dtsX   = zeros(Nship,1); dtPx=zeros(Nship,1);%
dtThe  = zeros(Nship,1); dtPthe =zeros(Nship,1);
dtExt  = zeros(Nship*6,1);

dz         = [dteta_hat; dtphi_hat;dtsX;dtPx; dtsZ; dtPz; dtThe; dtPthe; dtExt];


    if time>=par.t_init+(itersavevar-1)*par.dt
        if any (time<=shipsavevar.time(1:itersavevar-1))
            itersavevar=closest(shipsavevar.time,time);
        end
        
        shipsavevar.betaz(itersavevar,1:Nship)=beta.z(1:Nship);
        shipsavevar.betax(itersavevar,1:Nship)=beta.x(1:Nship);
        shipsavevar.phiwave(itersavevar,:)=phif;
        shipsavevar.time(itersavevar)=time;
        shipsavevar.zVel(itersavevar,1:Nship)=Vel.z;
        shipsavevar.xVel(itersavevar,1:Nship)=Vel.x;
        itersavevar=itersavevar+1;
    end


if par.IP.check==1
    if time>=par.IP.time(1) && time<=par.IP.time(2)
        if time>=par.IP.time(1)+(iterInteriordt-1)*par.dt*par.IP.time(3)
            dteta(iterInterior,:)=Ifft(dz(1:n));
            dtphihat(iterInterior,:)=dz(n+1:end)./(1i.*k);
            timeIP(iterInterior,1)=time;
            uIP(iterInterior,:)   =phi;
            etaIP(iterInterior,:) =eta;
            if any(k==0)
                dtphihat(isnan(dtphihat(iterInterior,:)))=0;
            end
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end

if strcmp(plot_ode,'yes')
    %if time>ITERdt*par.dt
        if ITERdt==1, fig2=figure;end;
        plot(x,eta,'r',x,phif,'--b',x,phitot,'-.g')
        title(['t=',num2str(time)]);
        ylim([-0.1 0.1])
        ITERdt=ITERdt+1;
        pause(0.001);
        
   % end
end
