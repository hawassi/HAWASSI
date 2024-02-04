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
%%%%%%%%%% RHS2D for varying bottom                                      %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dz] = RHSBathy2duv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,...
                           Oprt,timeSimul,ProgStatbar)

plotchecker=0;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=dom.fbl.char; %dampchar=char./length of fourier boundary;
dampcoef=7*Oprt.Cpeak;
aal     =dom.aal;
ChiAdj =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar  iterdtcheck Idstop
if ivp.type~=1
    if flagWarnIVP==1
        w = warning('query','last');  %turn off warning odeplot
        warning('off',w.identifier);
        flagWarnIVP=0;
    end
end

if par.IP.check==1
    global dteta dtphi iterInterior iterInteriordt timeIP phiIP etaIP
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
z_hat=z_hat.';
z1_hat=z_hat(1:N);      eta_hat=reshape(z1_hat,dom.Ny,[]);
z2_hat=z_hat(N+1:2*N);  u_hat=reshape(z2_hat,dom.Ny,[]);
z3_hat=z_hat(2*N+1:3*N);v_hat=reshape(z3_hat,dom.Ny,[]);

eta=funC_ifft2(eta_hat);u=funC_ifft2(u_hat);v=funC_ifft2(v_hat);
gradphi.x=u;gradphi.y=v;
gradphi_hat.x_hat=u_hat;gradphi_hat.y_hat=v_hat;
   
Depth  =-dom.bathy.profile ;   
H=(Depth+eta);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model.evol,'SWE')
if model.nonlinear==1  % still in phi formulation
depthgradphi_hat.x_hat=Depth.*gradphi_hat.x_hat;
depthgradphi_hat.y_hat=Depth.*gradphi_hat.y_hat;
div_depthgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,depthgradphi_hat);
delphiH_hat  = -div_depthgradphi_hat;
deletaH_hat  = g.*eta_hat;
else
Hgradphi_hat.x_hat=fft2(H.*gradphi.x);
Hgradphi_hat.y_hat=fft2(H.*gradphi.y);
div_Hgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,Hgradphi_hat);
gradphi2      = funOprt_innerproduct(gradphi,gradphi);
delphiH_hat  = -div_Hgradphi_hat;
deletaH_hat  = g.*eta_hat+0.5.*fft2(gradphi2);
end
dx_deletaH_hat =1i.*dom.Kx.*deletaH_hat;
dy_deletaH_hat =1i.*dom.Ky.*deletaH_hat;

elseif strcmp(model.evol,'HS')
Muv_hat=funOprt_M2d(g,dom.Kx,dom.Ky,Oprt.Cp2dSqInterp,u_hat,u,v_hat,v);
gradeta_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,eta_hat);    

delphiH1_hat     =Muv_hat;  % delphiH=-div deluH
dx_deletaH1_hat  =g.*gradeta_hat.x_hat;
dy_deletaH1_hat  =g.*gradeta_hat.y_hat;

if model.nonlinear==1
delphiH2_hat=0;    dx_deletaH2_hat=0; dy_deletaH2_hat=0;
delphiH3_hat=0;    dx_deletaH3_hat=0; dy_deletaH3_hat=0;   
else
M0uv                 = funC_ifft2(Muv_hat);
etaM0uv              = eta.*M0uv; etaMuv_hat= fft2(etaM0uv);
LetaM0uv_hat         = funOprt_L2d(g,Oprt.Om2dSqInterp,etaMuv_hat,etaM0uv);  
etagradphi_hat.x_hat=fft2(eta.*gradphi.x);
etagradphi_hat.y_hat=fft2(eta.*gradphi.y);
divetagradphi_hat   = funOprt_div2d(dom.Kx,dom.Ky,etagradphi_hat);
absuv2              = funOprt_innerproduct(gradphi,gradphi);
delphiH2_hat        = -(LetaM0uv_hat+divetagradphi_hat);
dx_deletaH2_hat     = 0.5.*1i.*dom.Kx.*fft2(absuv2-M0uv.^2);
dy_deletaH2_hat     = 0.5.*1i.*dom.Ky.*fft2(absuv2-M0uv.^2);

delphiH3_hat  =0; dx_deletaH3_hat=0;dy_deletaH3_hat=0;
  if model.nonlinear~=2
      LetaM0uv        =funC_ifft2(LetaM0uv_hat);
      etaLetaM0uv     =eta.*LetaM0uv;
      etaLetaM0uv_hat =fft2(etaLetaM0uv);
      LetaLetaM0uv_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,etaLetaM0uv_hat,etaLetaM0uv);
      eta2M0uv_hat    =fft2(eta.^2.*M0uv);
      gradeta2M0uv_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,eta2M0uv_hat);
      Divgradeta2M0uv_hat=funOprt_div2d(dom.Kx,dom.Ky,gradeta2M0uv_hat);
      DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,gradphi_hat);
      DivGradphi     =funC_ifft2(DivGradphi_hat);
      eta2Divgradphi =eta.^2.*DivGradphi;
      eta2Divgradphi_hat=fft2(eta2Divgradphi);
      Leta2DivGradphi_hat=funOprt_L2d(g,Oprt.Om2dSqInterp,eta2Divgradphi_hat,eta2Divgradphi);
      delphiH3_hat  =LetaLetaM0uv_hat+0.5.*(Divgradeta2M0uv_hat+Leta2DivGradphi_hat);
      dx_deletaH3_hat =1i.*dom.Kx.*fft2(M0uv.*(eta.*DivGradphi+LetaM0uv));
      dy_deletaH3_hat =1i.*dom.Ky.*fft2(M0uv.*(eta.*DivGradphi+LetaM0uv));
  end
end
delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat+delphiH3_hat).*ChiAdj);
dx_deletaH_hat=dx_deletaH1_hat+fft2(funC_ifft2(dx_deletaH2_hat+dx_deletaH3_hat).*ChiAdj);
dy_deletaH_hat=dy_deletaH1_hat+fft2(funC_ifft2(dy_deletaH2_hat+dy_deletaH3_hat).*ChiAdj);
end

if strcmp(model.breaking.check,'Yes')
    [B] = breaking_process2d(eta,gradphi,dom,parBreak,time,g,t_init);
    visc=B.*parBreak.delb.^2.*H.*funC_ifft2(delphiH_hat);
    Flux=-visc.*funC_ifft2(delphiH_hat);
    [dxFlux, dyFlux]=gradient(Flux,dom.dx,dom.dy);
    Sb_hat.x  = fft2(dxFlux./H);
    Sb_hat.y  = fft2(dyFlux./H);    
else
    Sb_hat.x=0; Sb_hat.y=0;
end

if dom.friction.input.check==1
    absGradphi=sqrt(gradphi.x.^2+gradphi.y.^2);
    Sfx=-dom.friction.Char.*(1./(H.^(4/3))).*u.*absGradphi;
    Sfy=-dom.friction.Char.*(1./(H.^(4/3))).*v.*absGradphi;
    Sf_hat.x=fft2(Sfx);
    Sf_hat.y=fft2(Sfy);
else, Sf_hat.x=0;Sf_hat.y=0;
end

dteta_hat=(delphiH_hat +Source_hat-dampcoef*fft2(eta.*dampchar)).*aal;
dtu_hat=(-dx_deletaH_hat+Sb_hat.x+Sf_hat.x-dampcoef*fft2(u.*dampchar)).*aal;
dtv_hat=(-dy_deletaH_hat+Sb_hat.y+Sf_hat.y-dampcoef*fft2(v.*dampchar)).*aal;

dz1_hat=reshape(dteta_hat,1,[]);
dz2_hat=reshape(dtu_hat,1,[]);
dz3_hat=reshape(dtv_hat,1,[]);
dz=[dz1_hat,dz2_hat,dz3_hat].';


if par.IP.check==1
    if time>par.IP.time(1) && time<par.IP.time(2)
        if time>par.IP.time(1)+(iterInteriordt-1)*dt*par.IP.time(3)
            dteta(iterInterior,:,:)   =funC_ifft2(dteta_hat);
            dtphi(iterInterior,:,:)   =intgrad2(funC_ifft2(dtu_hat),funC_ifft2(dtv_hat),dom.dx,dom.dy);
            timeIP(iterInterior,1)    =time;
            phiIP(iterInterior,:,:)   =intgrad2(u,v,dom.dx,dom.dy);
            etaIP(iterInterior,:,:)   =eta;
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end


if plotchecker==1
[X,Y]=meshgrid(dom.X,dom.Y);
if  iterdtcheck==1
    figure;
end
max(max(eta))
if mod( iterdtcheck,2)
surf(X,Y,eta,'Edgecolor','none')
view(2)
colorbar;
caxis([-0.1;0.1]);
zlim([-0.1;0.1])
title(['time= ',num2str(time)])
pause(0.001);
end
 iterdtcheck= iterdtcheck+1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if max(any(abs(eta)>6*par.Hs))==1 || min(min(eta))< min(min(-Depth))%%
%     Idstop=2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
