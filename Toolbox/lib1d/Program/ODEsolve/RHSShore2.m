% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%    HAWASSI 1D                                          %%%%%%%%%%
%%%%%%%%%    Hamiltonian Wave-Ship-Structure Interaction         %%%%%%%%%%
%%%%%%%%%    copyright (c): LabMath-Indonesia & Unversity Twente %%%%%%%%%%
%%%%%%%%%    version 23 June 2015                                %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  RHS for RunUp  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dz = RHSShore2(time,z,model,par,influx,IVP,Oprt,ProgStatbar)
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
H_min       =par.interp.H_min;
H_mid       =par.interp.H_mid;
H_max       =par.interp.H_max;
IntCoef     =par.interp.IntCoef;
Up          =Oprt.fun_Up;
omAdd       = model.OmFun;
Cder        =Oprt.fun_Cder;

  

global flagWarnIVP iterProgBar ITERdt
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
end;


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

if isfield(par,'wall')
if strcmp(par.wall.presence,'Yes')
Wallchar    =par.wall.Rho';
indEndBreak =closest(x,par.wall.position);
else
Wallchar    =ones(n,1);  
indEndBreak =n;
end
end 

H=Chi.*eta+depth;
Upeta   = Ifft(Oprt.Cpeak.*z1);
Upu     = Ifft(Oprt.Cpeak.*z2);

HeavH=Heaviside(H-H_min);

indW=find(HeavH==1);
indD=find(HeavH<H_min);
ShoreChar=HeavH;  



if strcmp(model.dispersion,'OmSWE')   % Use direct model to avoid interpolation
    if model.nonlinear==2
    Htot        =(depth+Chi.*eta);
    deluH_hat   =fft(Htot.*Ifft(z2).*ShoreChar); 
    deletaH_hat = fft(Ifft(g*z1+fft(Chi.*u.^2)./2).*ShoreChar);
    deta_hat     = -1i*k.*deluH_hat ;
    du_hat       = -1i*k.*deletaH_hat ;
    end
elseif strcmp(model.evol,'HS')    
    Cm         = Up(k,H_min,omAdd); ConjCm=conj(Cm);
    Cp         = Up(k,H_max,omAdd); ConjCp=conj(Cp);
    Cc         = Up(k,H_mid,omAdd); ConjCc=conj(Cc);
    Cder_m     = Cder(k,H_min);
    Cder_p     = Cder(k,H_max);
    Cder_c     = Cder(k,H_mid);
        
    gm=zeros(n,1);gp=gm;gc=gp;

    gm(indW)=ppval(IntCoef.m,H(indW));
    gp(indW)=ppval(IntCoef.p,H(indW));
    gc(indW)=ppval(IntCoef.c,H(indW));

    
    Cu_hat                = Op_C(Cm,gm,Cp,gp,Cc,gc,u,z2);
    Cu                    = Ifft(Cu_hat);
    Cu_ShoreChar          = Ifft(Cu_hat).*ShoreChar.*Wallchar;
    Cu_ShoreChar_hat      = fft(Cu_ShoreChar);
    ConjC_Cu_ShoreChar_hat= Op_C(ConjCm,gm,ConjCp,gp,ConjCc,gc,...
    Cu_ShoreChar,Cu_ShoreChar_hat);
    Cder_u_hat             = Op_C(Cder_m,gm,Cder_p,gp,Cder_c,gc,u,z2);
    Cder_u                = Ifft(Cder_u_hat);
    
    deluH_hat   =  1/g.*ConjC_Cu_ShoreChar_hat;
    deletaH_hat =  g.*z1+ fft(Ifft(fft(Cder_u.*Cu)./g).*Chi.*Wallchar);
    deletaH_hat = fft(Ifft(deletaH_hat).*ShoreChar);
    
    deta_hat = -1i.*k.*(deluH_hat);
    du_hat   = -1i.*k.*(deletaH_hat);    
end




%Breaking process
if strcmp(model.breaking.check,'Yes')  
    [B]=breaking_process(eta,u,depth,time,x,par.break,indEndBreak);
    visc=B.*par.break.delb.^2.*(H).*ifft(deta_hat,'symmetric');
    Flux=-visc.*Ifft(deta_hat);
    dxFlux=gradient(Flux,par.dx);
    oneperH=zeros(n,1);
    oneperH(H>H_min)=1./H(H>H_min);
    Rb_hat=fft(oneperH.*dxFlux.*Wallchar);
    
        if time>28
        subplot(2,1,2)
        Rb=Ifft(Rb_hat);
        plot(x,Rb,'c',x,visc,'--r',x,oneperH,':b',x,u,'k',x,ShoreChar,'g--')
        title(['time= ', num2str(time)])
        xlim([20;100]);ylim([-10;10]);
        pause(0.001)
        end
else
    Rb_hat=0;
end

if par.friction.check==1
Cm=par.friction.Cf;%%manning coefficient [0.01;0.04]
indexf=par.friction.index;
Rf=zeros(n,1);
Rf(indexf)=-g*Cm.^2.*(u(indexf).*abs(u(indexf))./(H(indexf)).^(4/3));
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


signProp=1;

dampchar(indD)=ones;
Upeta(indD)   =H(indD);
Upu (indD)    =u(indD);


deta_hat = (deta_hat + Source -(1+signProp*1i.*k).*fft(Upeta.*dampchar)).*aal;
du_hat   = (du_hat+Rb_hat+Rf_hat+Rw_hat-(1+signProp*1i.*k).*fft(Upu.*dampchar)).*aal;

dz  = [deta_hat ; du_hat];

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

     

   


if strcmp(Visual_ODE,'On')
    if time>ITERdt*par.dt
        %if ITERdt==1, fig2=figure;end;
        plot(x,eta,'r',x,par.bathy,'--k',x,ShoreChar,'b:')
        title(['t=',num2str(time)]);
       % ylim([-0.2 0.2])
        ITERdt=ITERdt+1;
        pause(0.0001)
        
    end
end

