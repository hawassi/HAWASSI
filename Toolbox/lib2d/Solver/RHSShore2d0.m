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
%%%%%%%%%% RHS2D for flat bottom                                      %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dz] = RHSShore2d(time,z_hat,input,model,par,parBreak,dom,influx,ivp,...
                timeSimul,Oprt,ProgStatbar)

plotchecker=1;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=dom.fbl.char; %dampchar=char./length of fourier boundary;
dampcoef=7*Oprt.Cpeak;
aal     =dom.aal;
Depth  =-dom.bathy.profile ;
ChiAdj =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
H_min=Oprt.H_min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar  iterdtcheck
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
z_hat=z_hat.';
z1_hat=z_hat(1:N);     eta_hat=reshape(z1_hat,dom.Ny,[]);
z2_hat=z_hat(N+1:end); phi_hat=reshape(z2_hat,dom.Ny,[]);
eta=funC_ifft2(eta_hat);phi=funC_ifft2(phi_hat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Define active domain by a characteristic function
WaveChar=ones(size(eta));
H=ChiAdj.*eta+Depth;
WaveChar(H<H_min)=zeros;
IndxZerofirst=funC_findIndexinColumn(WaveChar,'=',0,'first');
WaveChar(:,IndxZerofirst(:):end)=0;
dampchar(WaveChar==0)=ones;
%dampchar(:,IndxZerofirst(:):end)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model.evol,'SWE')

% if model.nonlinear==1
% gradphi_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,phi_hat);    
% depthgradphi_hat.x_hat=Depth.*gradphi_hat.x_hat;
% depthgradphi_hat.y_hat=Depth.*gradphi_hat.y_hat;
% div_depthgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,depthgradphi_hat);
% delphiH_hat  = -div_depthgradphi_hat;
% deletaH_hat  = g.*eta_hat;
% else
% gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);      
% Hgradphi_hat.x_hat=fft2((Depth+eta).*gradphi.x);
% Hgradphi_hat.y_hat=fft2((Depth+eta).*gradphi.y);
% div_Hgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,Hgradphi_hat);
% gradphi2      = funOprt_innerproduct(gradphi,gradphi);
% delphiH_hat  = -div_Hgradphi_hat;
% deletaH_hat  = g.*eta_hat+0.5.*fft2(gradphi2);
% end;


elseif strcmp(model.evol,'HS')
gradphi              = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
gradphixChar         =gradphi.x.*WaveChar;gradphiyChar=gradphi.y.*WaveChar;

Cp2gradphiChar_hat.x_hat = funOprt_Cp2d(Oprt.Cp2dSqInterp,...
                             fft2(gradphixChar),gradphixChar);
Cp2gradphiChar.x         =   funC_ifft2(Cp2gradphiChar_hat.x_hat);
Cp2gradphiChar.x(Depth<0)=-1.*Cp2gradphiChar.x(Depth<0);
Cp2gradphiChar_hat.x_hat =fft2(Cp2gradphiChar.x);

Cp2gradphiChar_hat.y_hat = funOprt_Cp2d(Oprt.Cp2dSqInterp,...
                             fft2(gradphiyChar),gradphiyChar);
Cp2gradphiChar.y         =   funC_ifft2(Cp2gradphiChar_hat.y_hat);
Cp2gradphiChar.y(Depth<0)=-Cp2gradphiChar.y(Depth<0);
Cp2gradphiChar_hat.y_hat =fft2(Cp2gradphiChar.y);

LphiChar_hat=-funOprt_div2d(dom.Kx,dom.Ky,Cp2gradphiChar_hat)/g;

delphiH1_hat  = LphiChar_hat;
deletaH1_hat  =g.*fft2(eta.*WaveChar);
if model.nonlinear==2
Lphi              = funC_ifft2(funOprt_L2d(Oprt.Om2dSqInterp,phi_hat,phi));
Lphi(Depth<0)     =-Lphi(Depth<0);
etaLphiChar       = eta.*Lphi.*WaveChar;
etaLphiChar_hat   = fft2(etaLphiChar);
LetaLphiChar_hat  = funOprt_L2d(Oprt.Om2dSqInterp,etaLphiChar_hat,etaLphiChar); 
LetaLphiChar      =funC_ifft2(LetaLphiChar_hat );
LetaLphiChar(Depth<0)=-LetaLphiChar(Depth<0);
LetaLphiChar_hat  =fft2(LetaLphiChar);
etagradphiChar_hat.x_hat=fft2(eta.*gradphixChar);
etagradphiChar_hat.y_hat=fft2(eta.*gradphiyChar);
divetagradphiChar_hat = funOprt_div2d(dom.Kx,dom.Ky,etagradphiChar_hat);
gradphi2      = funOprt_innerproduct(gradphi,gradphi);
delphiH2_hat  = -(LetaLphiChar_hat+divetagradphiChar_hat);
deletaH2_hat  = 0.5*fft2((gradphi2-Lphi.^2).*WaveChar); 
else
delphiH2_hat=0;    deletaH2_hat=0;
end
delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat).*ChiAdj);
deletaH_hat=deletaH1_hat+fft2(funC_ifft2(deletaH2_hat).*ChiAdj);
end

if strcmp(model.breaking.check,'Yes')
    [B] = breaking_process2d(eta,gradphi,dom,parBreak,time,t_init);
    Sb_hat  = -parBreak.delb.^2*fft2(B.*funC_ifft2(delphiH_hat).^2);
else
    Sb_hat=0;
end

if dom.friction.input.check==1
    Sf=-dom.friction.Char.*(1./H).*phi;
    Sf(H<H_min)=0;Sf_hat=fft2(Sf);
else Sf_hat=0;
end
etaH=eta;
etaH(WaveChar==0)=H(WaveChar==0);
dteta_hat=(delphiH_hat+Source_hat-dampcoef.*(1+1i.*dom.KK).*fft2(etaH.*dampchar)).*aal;
dtphi_hat=(-deletaH_hat+Sb_hat+Sf_hat-dampcoef.*(1+1i.*dom.KK).*fft2(phi.*dampchar)).*aal;

dz1_hat=reshape(dteta_hat,1,[]);
dz2_hat=reshape(dtphi_hat,1,[]);
dz=[dz1_hat,dz2_hat].';


if plotchecker==1
stepS=4;
[X,Y]=meshgrid(dom.X(1:stepS:end),dom.Y(1:stepS:end));
if  iterdtcheck==1
    figure;
end
if mod( iterdtcheck,2)
hold off;
surf(X,Y,eta(1:stepS:end,1:stepS:end),'Edgecolor','none')
% hold on
% surf(X,Y,WaveChar(1:stepS:end,1:stepS:end),'Edgecolor','none')
view([-10,10])
colorbar;
caxis([-0.05;0.05]);
zlim([-0.05;0.05])
title(['time= ',num2str(time)])
pause(0.01);
end
 iterdtcheck= iterdtcheck+1;
end
