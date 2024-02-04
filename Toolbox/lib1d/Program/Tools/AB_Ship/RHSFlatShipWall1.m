%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 12 Dec 2019                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for Flat-Bottom %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dz] = RHSFlatShipWall1(time,z,model,par,influx,IVP,Oprt,shippar,ProgStatbar)
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
aal          =Oprt.aal;
aalS         =Oprt.aalS;
L0           =Oprt.L0;
C0           =Oprt.Upd;
Wallchar    =par.wall.Rho';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar ITERdt shipsavevar itersavevar Sxprev tprev
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
sTheta  = z(2*n+4*Nship+1:2*n+Nship*5);
Ptheta  = z(2*n+5*Nship+1:2*n+Nship*6);
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
kappaL      =shippar.difpot.kappaL; 
kappaR      =shippar.difpot.kappaR;

Amat       =shippar.difpot.Amat;
Bvect      =shippar.difpot.Bvect;
Smotion    =shippar.data(:,2);

%     indl=closest(x,shippar.form.xShip(1,1));indr=closest(x,shippar.form.xShip(1,3));
%     Xcnow=shippar.form.xShip(1,2);
%     IdMove=0;
%     if Xcnow<Sxprev(1,2)
% %         Nshiftl=indl-indlprev;
% %         phi(indl:indl-Nshiftl)=phi(indlprev:indlprev-Nshiftl);
%             
%         phiLnow=phi(indl);phiLprev=phiprev(indlprev);phiRprev=phiprev(indrprev);
%         phi(indr)=phiRprev+(phiLnow-phiLprev);
%         phi(indr+1:indrprev-1)=phi(indr)+(x(indr+1:indrprev-1)-x(indr)).*(phi(indrprev)-phi(indr))./(x(indrprev)-x(indr));
%         IdMove=-1;
%     elseif Xcnow>Sxprev(1,2)
% %         Nshiftr=indr-indrprev;
% %         phi(indr-Nshiftr:indr)=phi(indrprev-Nshiftr:indrprev);
%         
%         phiRnow=phi(indr);phiRprev=phiprev(indrprev);phiLprev=phiprev(indlprev);
%         phi(indl)=phiLprev+(phiRnow-phiRprev);
%         phi(indlprev+1:indl-1)=phi(indl)+(x(indlprev+1:indl-1)-x(indl)).*(phi(indlprev)-phi(indl))./(x(indlprev)-x(indl));
%         IdMove=1;
%     end

if strcmp(model.evol,'HS')
   if shippar.calclin==1
       sZInDif=0;
   else
       sZInDif=sZ;
   end
        [phif,Fs,FEv,Lev,beta,dxi_Kdiff,etaL,etaR,dxetaL,dxetaR]=funS_tangential_potential3...
        (shippar.calclin,k,x,xShip,Sdraft,Slength,sZInDif,phi,eta,chiship,chiwaterline,...
        -par.bathy',Nship,Evmode,kappaL,kappaR,Amat,Bvect,nutild,shippar.form.XZrot0(1,:));
   
    if shippar.calclin==0
        if any(~strcmp(Smotion(:),'Fixed'))
            [shippar.rad.psi,shippar.rad.Ma,shippar.rad.dampC,shippar.rad.LpsiEv,...
                shippar.rad.Fpsi,shippar.rad.FpsiEv,shippar.rad.FpsiS,dxi_Krad]=...
                funS_radiation_potential(shippar,x,xShip,Sdraft,sZ,chiship,...
                chiwaterline,nu,Smotion,-par.bathy',Nship,etaL*0,etaR*0,dxetaL*0,dxetaR*0);        
        else
            dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
         end
    else
        dxi_Krad.z(1)=0;dxi_Krad.x(1)=0;dxi_Krad.theta(1)=0;
    end
    [phiRad,Vel,dtetaS,FpsiS,dxi_Krad]=fun_vel_and_radpot_calc(shippar,chiship,...
         nutild,Pz,Px,Ptheta,beta,Smotion,Nship,dxi_Krad);
    
    phitot=phif+phiRad;
    dphidtS=((phitot-phiprev)/par.dt);
      
    FluxF  =ifft(funOprt_Flux(g,k,C0,eta,phi,model.nonlinear),'symmetric');
    
    Lphif_hat=funOprt_LWall(g,k,C0,eta,phi,model.nonlinear,ChiAdj,Wallchar);
    
    chishipF=chiship(:,end);
   
    deletaH_hat  = fft(g*eta);
    if model.nonlinear~=1
        phitot_hat=fft(phi);
        dxphi=Ifft(1i.*k.*phitot_hat);
        L0phi=Ifft(L0.*phitot_hat);
        deletaH_hat  = (deletaH_hat+0.5.*fft((dxphi.^2-L0phi.^2).*ChiAdj));
        if model.nonlinear~=2
            deletaH_hat  = (deletaH_hat+ fft((L0phi.*(eta.*Ifft(1i.*k.*fft(dxphi))+...
                Ifft(L0.*fft(eta.*L0phi)))).*ChiAdj));
        end
    end
    
%     if IdMove==1 || IdMove==-1
%     aal=aalS;
%     chishipPrev=zeros(size(eta));
%     chishipPrev(indlprev+1:indrprev-1)=1;
%     chishipDelMove=chishipF-chishipPrev;
%      etaN=eta;
%     if IdMove==1
%         chishipDelMove(indrprev)=0;
%          etaN(indl)=eta(indl-1);
%      elseif IdMove==-1
%         chishipDelMove(indlprev)=0; 
%         etaN(indr)=eta(indr+1);
%      end
%      detadtS=(etaN-etaprev)./par.dt;
%     else
%       chishipDelMove=0;
%     end
    LphiF_hat=fft(Ifft((Lphif_hat+fft(Lev)+fft((FluxF+FEv+-Fs-FpsiS).*dxchiship(:,end))).*aal).*chifluid+dtetaS);
    dteta_hat     =   LphiF_hat; 
    dtphi_hat     = -fft(Ifft(deletaH_hat).*chifluid-dphidtS.*chishipF).*aal ;
    
end

if strcmp(model.breaking.check,'Yes')
     u=Ifft(1i.*k.*fft(phi));
%     [B] = breaking_process(eta,u,depth,time,x,par.break,[]);
%     Rb_hat  = -par.break.delb.^2.*fft(B.*Ifft(dteta_hat).^2);
[B] = breaking_process(eta,u,depth,time,x,par.break,[]);
    visc= B.*par.break.delb.^2.*(depth+eta).*Ifft(dteta_hat);
    Rb_hat  = -fft((1./(depth+eta)).*Ifft(1i.*k.* fft(visc.*Ifft(dteta_hat)) ))./(1i.*k);
 
else
    Rb_hat=0;
end

if par.friction.check==1
    H=depth+eta;
    Cf=par.friction.Cf;%%chezy coefficient [0.001;0.01]
    indexf=par.friction.index;
    u=Ifft(1i.*k.*fft(phi));
    Rf=zeros(n,1);
    Rf(indexf)=-Cf.*(u(indexf).*abs(u(indexf))./H(indexf));
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

dteta_hat   = dteta_hat + (Source   - fft(par.dampcoef.*eta.*dampchar)).*aal;
dtphi_hat     = dtphi_hat +(Rb_hat +Rf_hat+Rw_hat  - fft(par.dampcoef.*phi.*dampchar)).*aal;


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
    dtPimoored=((shippar.form.MassMat+0*shippar.rad.Ma.Mat)*xi).*(2*pi./Tn).^2;
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
        
        if strcmpi(shippar.moored.check,'yes')&& Tn(1)~=0
            dtPx(ii,1)   = dtPx(ii,1)-dtPimoored(1)-b.*xi(1);
        end
%         if shippar.calclin==0
%             dtPx(ii,1)   = dtPx(ii,1) -dxi_Kdiff.x(1) -dxi_Krad.x(1);
%         end
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
%         if shippar.calclin==0
%             dtPz(ii,1)   =dtPz(ii,1)  -dxi_Kdiff.z(1)-dxi_Krad.z(1);
%         end
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
        if strcmpi(shippar.moored.check,'yes')&& Tn(3)~=0
            dtPthe(ii,1)   = dtPthe(ii,1)-dtPimoored(3) ;
        end
        
%         if shippar.calclin==0
%             dtPthe(ii,1)   =dtPthe(ii,1) -dxi_Kdiff.theta(1)-dxi_Krad.theta(1);
%         end
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
    shipsavevar.addMass(itersavevar,:)=reshape(shippar.rad.Ma.Mat,1,[]);
    shipsavevar.dampCoef(itersavevar,:)=reshape(shippar.rad.dampC,1,[]);
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
     if time>ITERdt*par.dt && time>0
        if ITERdt==1, fig2=figure;end
        
%         chishipOld    =chiship;
%         xS=sX(1)+shippar.form.xShip(1,2);
%         indlprev=closest(x,Sxprev(1,1));indrprev=closest(x,Sxprev(1,3));
%         [shippar.form]= shapeformsetup(x,shippar,xS,sZ(1));
%         chishipNew    =shippar.form.chi;
%         
%          
%         phiprev=phi;  
%        % eta=eta.*chifluid;
%         indl=closest(x,shippar.form.xShip(1,1));indr=closest(x,shippar.form.xShip(1,3));
%         Xcnow=shippar.form.xShip(1,2);
%         
%         %phi(indl:indr)=phiprev(indlprev:indrprev);
%             
%         if Xcnow<Sxprev(1,2)
%           %  Nshiftl=indl-indlprev
%             %phi(indl:indl-Nshiftl)=phi(indlprev:indlprev-Nshiftl);
%             
%             phiLnow=phi(indl);phiLprev=phiprev(indlprev);phiRprev=phiprev(indrprev);
%             phi(indr)=phiRprev+(phiLnow-phiLprev);
%             phi(indr+1:indrprev-1)=phi(indr)+(x(indr+1:indrprev-1)-x(indr)).*(phi(indrprev)-phi(indr))./(x(indrprev)-x(indr));
% %           phi=smooth(phi,2);
%         elseif Xcnow>Sxprev(1,2)
% %             Nshiftr=indr-indrprev
% %             phi(indr-Nshiftr:indr)=phi(indrprev-Nshiftr:indrprev);
% %             
%             phiRnow=phi(indr);phiRprev=phiprev(indrprev);phiLprev=phiprev(indlprev);
%             phi(indl)=phiLprev+(phiRnow-phiRprev);
%             phi(indlprev+1:indl-1)=phi(indl)+(x(indlprev+1:indl-1)-x(indl)).*(phi(indlprev)-phi(indl))./(x(indlprev)-x(indl));
% %             phi=smooth(phi,2);
%         end
%         
% 
% plot(x,chishipOld,'b',x,chishipNew,'r',x,phiprev,'b--',x,phi,'r--');
% title(['Xprev= ',num2str(Sxprev(1,2)),' Xnow=',num2str(Xcnow)])
% xlim([-60;60])
% pause;
plot(x,eta)
       pause(0.001)  
        ITERdt=ITERdt+1;
         
    end
end
 Sxprev=shippar.form.xShip;

