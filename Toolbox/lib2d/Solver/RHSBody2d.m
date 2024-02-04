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
%%%%%%%%%% HAWASSI-AB 2D                                         %%%%%%%%%%
%%%%%%%%%% Hamiltonian Wave-Ship-Structure Interaction           %%%%%%%%%%
%%%%%%%%%% Copyright (c): LabMath-Indonesia                      %%%%%%%%%%
%%%%%%%%%% version: 5 July 2016                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% RHS2D for wave-floating body interaction              %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dz] = RHSBody2d(time,z,input,model,par,dom,influx,ivp,body,Oprt,...
                          timeSimul,ProgStatbar)

plotchecker=0;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=dom.fbl.char; %dampchar=char./length of fourier boundary;
dampcoef=7*Oprt.Cpeak;
aal     =dom.aal;
ChiAdj =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
%%%%%%%%parameters of a body%%%%%%%%%%%%%%%%%%%%%%%%
Nbody=body.N;
Mass =body.par.mass;
aMz  =body.addedmass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar  iterdtcheck bodysavevar itersaveBodyvar
if ivp.type~=1
    if flagWarnIVP==1
        w = warning('query','last');  %turn off warning odeplot
        warning('off',w.identifier);
        flagWarnIVP=0;
    end
end

%%Progress bar
t_init=timeSimul.t_init; dt=timeSimul.dt;
Nt=length(timeSimul.interval);
if time>t_init+iterProgBar*dt*floor(0.1*Nt) %update every 10%
    set(jProgressBar,'Maximum',Nt, 'Value',iterProgBar*floor(0.1*Nt));
    jProgressBar.setStringPainted( true );
    ETA=funGui_remain_time(iterProgBar*floor(0.1*Nt),Nt);
    statusbarObj.setText(['time remaining=', num2str(ETA)]);
    iterProgBar=iterProgBar+1;
end
%reshape matrix Ny.Nx
N=dom.Nx*dom.Ny;
z=z.';
z1_hat=z(1:N);     eta_hat=reshape(z1_hat,dom.Ny,[]);
z2_hat=z(N+1:2*N); phi_hat=reshape(z2_hat,dom.Ny,[]);
eta=funC_ifft2(eta_hat);phi=funC_ifft2(phi_hat);

% sX          = z(2*N+1:2*N+Nbody);
% sY          = z(2*N+Nbody+1:2*N+Nbody*2);
% sZ          = z(2*N+2*Nbody+1:2*N+Nbody*3);
% Degx        = z(2*N+3*Nbody+1:2*N+Nbody*4);
% Degy        = z(2*N+4*Nbody+1:2*N+Nbody*5);
% Degz        = z(2*N+5*Nbody+1:2*N+Nbody*6);
% PsX         = z(2*N+6*Nbody+1:2*N+Nbody*7);
% PsY         = z(2*N+7*Nbody+1:2*N+Nbody*8);
 PsZ         = z(2*N+8*Nbody+1:2*N+Nbody*9);
% PDegX       = z(2*N+9*Nbody+1:2*N+Nbody*10);
% PDegY       = z(2*N+10*Nbody+1:2*N+Nbody*11);
% PDegZ       = z(2*N+1*Nbody+1:2*N+Nbody*12);
% xVel        =zeros(Nbody,1);
% yVel        =zeros(Nbody,1);
zVel        =zeros(Nbody,1);
betaz       =zeros(Nbody,1);
betax       =zeros(Nbody,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%parameters%%%%%%%%%%
chiB=body.par.chiB(end).char;chiF=1-chiB;

zeta=eta.*chiF+body.par.shape.*chiB;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%phitan=TangentialPotential2dbyMinimization(dom,Oprt,chiB,phi,phi_hat);
phitan=TangentialPotential2dbyMinimization1(g,dom,Oprt,chiB,phi,phi_hat,zeta);
%phitan=TangentialPotential2d(dom,Oprt,chiB,phi,phi_hat);

phirad=0;LphiS=0;ZetaSI(1:Nbody).profile=zeros;
for i=1:Nbody
chiBI=body.par.chiB(i).char;
betaz(i)    = trapz(dom.Y,trapz(dom.X,phitan.*chiBI,2));
betax(i)    = trapz(dom.Y,trapz(dom.X,phitan.*body.par.dxshape.*chiBI,2));
if body.fixed==0
zVel(i)     = (PsZ(i)-betaz(i))/(Mass(i)+aMz(i));
else
zVel(i)=0;    
end
phiRad=phirad+zVel(i)*body.radpot.*chiBI;
LphiS=LphiS+zVel(i).*chiBI;
ZetaSI(i).profile=zeta.*chiBI;
end

phitot=phitan+phiRad;

if strcmp(model.evol,'SWE')

elseif strcmp(model.evol,'HS')
Lphi_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,fft2(phitot),phitot);
LphiFS_hat=fft2(funC_ifft2(Lphi_hat).*chiF+LphiS);
dphidtS=(phitot-phi)/(timeSimul.dt);
delphiH1_hat  = LphiFS_hat;
deletaH1_hat  =fft2(g.*eta.*chiF-dphidtS.*chiB);
if model.nonlinear==2
% Lphi          = funC_ifft2(Lphi_hat);
% etaLphi       = eta.*Lphi;
% etaLphi_hat   = fft2(etaLphi);
% LetaLphi_hat  = funOprt_L2d(Oprt.Om2dSqInterp,etaLphi_hat,etaLphi);  
% gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);   
% etagradphi_hat.x_hat=fft2(eta.*gradphi.x);
% etagradphi_hat.y_hat=fft2(eta.*gradphi.y);
% divetagradphi_hat = funOprt_div2d(dom.Kx,dom.Ky,etagradphi_hat);
% gradphi2      = funOprt_innerproduct(gradphi,gradphi);
% delphiH2_hat  = -(LetaLphi_hat+divetagradphi_hat);
% deletaH2_hat  = 0.5*fft2(gradphi2-Lphi.^2); 
else
delphiH2_hat=0;    deletaH2_hat=0;
end
delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat).*ChiAdj);
deletaH_hat=deletaH1_hat+fft2(funC_ifft2(deletaH2_hat).*ChiAdj);
end

dteta_hat=(delphiH_hat+Source_hat-dampcoef*fft2(eta.*dampchar)).*aal;
dtphi_hat=(-deletaH_hat-dampcoef*fft2(phi.*dampchar)).*aal;

dz1_hat=reshape(dteta_hat,1,[]);
dz2_hat=reshape(dtphi_hat,1,[]);
dSxyz=zeros(3*Nbody,1);
dSxyz(2*Nbody+1:3*Nbody)=zVel;
dDegxyz=zeros(3*Nbody,1);
dPSD=zeros(6*Nbody,1);
dPSD(2*Nbody+1:3*Nbody)=-g*Mass-g.*trapz(dom.Y,trapz(dom.X,ZetaSI(i).profile,2));

dz=[dz1_hat,dz2_hat,dSxyz',dDegxyz',dPSD'].';

if body.saveoutputvar==1
   if time>=influx.timesig(1)+(itersaveBodyvar-1)*influx.par.dt
        if any (time<=bodysavevar.time(1:itersaveBodyvar-1))
            itersaveBodyvar=funC_closest(bodysavevar.time,time);
        end
        bodysavevar.betaz(itersaveBodyvar,1:body.N)=betaz(1:body.N);
        bodysavevar.betax(itersaveBodyvar,1:body.N)=betax(1:body.N);
        bodysavevar.time(itersaveBodyvar)=time;
        itersaveBodyvar=itersaveBodyvar+1;
    end 
    
end

if plotchecker==1
StepS=4;
[X,Y]=meshgrid(dom.X(1:StepS:end),dom.Y(1:StepS:end));
if  iterdtcheck==1
    figure;
end

if mod( iterdtcheck,2)
%subplot(3,1,1)
surf(X,Y,zeta(1:StepS:end,1:StepS:end),'Edgecolor','none')
view([10,10])
colorbar;
caxis([-0.3;0.1]);
zlim([-0.3;0.1])
title(['time= ',num2str(time)])
% subplot(2,1,2)
% surf(X,Y,phitan,'Edgecolor','none')
% view(2)
% colorbar;
%caxis([-0.1;0.1]);
%zlim([-0.1;0.1])
pause(0.01);
end
 iterdtcheck= iterdtcheck+1;
end
