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
%%%%%%%%%% RHS2D for flat bottom                                 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dz] = RHSFlat2d(time,z_hat,input,model,par,parBreak,dom,influx,ivp,...
    Oprt,timeSimul,ProgStatbar)

plotchecker=0;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=dom.fbl.char; %dampchar=char./length of fourier boundary;
dampcoef=7*Oprt.Cpeak;
aal     =dom.aal;
Depth  =-dom.bathy.profile ;
ChiAdj =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;
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
    global dteta dtphihat iterInterior iterInteriordt timeIP phiIP etaIP
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model.evol,'SWE')
    
    if model.nonlinear==1
        gradphi_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,phi_hat);
        depthgradphi_hat.x_hat=Depth.*gradphi_hat.x_hat;
        depthgradphi_hat.y_hat=Depth.*gradphi_hat.y_hat;
        div_depthgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,depthgradphi_hat);
        delphiH_hat  = -div_depthgradphi_hat;
        deletaH_hat  = g.*eta_hat;
    else
        gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
        Hgradphi_hat.x_hat=fft2((Depth+eta).*gradphi.x);
        Hgradphi_hat.y_hat=fft2((Depth+eta).*gradphi.y);
        div_Hgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,Hgradphi_hat);
        gradphi2      = funOprt_innerproduct(gradphi,gradphi);
        delphiH_hat  = -div_Hgradphi_hat;
        deletaH_hat  = g.*eta_hat+0.5.*fft2(gradphi2);
    end
    
    
elseif strcmp(model.evol,'HS')
    Lphi_hat=Oprt.L2d.*phi_hat;
    
    delphiH1_hat  = Lphi_hat;
    deletaH1_hat  =g.*eta_hat;
    if model.nonlinear==1
        delphiH2_hat=0;    deletaH2_hat=0;
        delphiH3_hat=0;    deletaH3_hat=0;
    else
        Lphi          = funC_ifft2(Lphi_hat);
        LetaLphi_hat  = Oprt.L2d.*fft2(eta.*Lphi);
        gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
        etagradphi_hat.x_hat=fft2(eta.*gradphi.x);
        etagradphi_hat.y_hat=fft2(eta.*gradphi.y);
        divetagradphi_hat = funOprt_div2d(dom.Kx,dom.Ky,etagradphi_hat);
        gradphi2      = funOprt_innerproduct(gradphi,gradphi);
        delphiH2_hat  = -(LetaLphi_hat+divetagradphi_hat);
        deletaH2_hat  = 0.5*fft2(gradphi2-Lphi.^2);
        delphiH3_hat  =0; deletaH3_hat=0;
        if model.nonlinear~=2
            LetaLphi        =funC_ifft2(LetaLphi_hat);
            LetaLetaLphi_hat=Oprt.L2d.*fft2(eta.*LetaLphi);
            eta2Lphi_hat=fft2(eta.^2.*Lphi);
            gradeta2Lphi_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,eta2Lphi_hat);
            DivGradeta2Lphi_hat=funOprt_div2d(dom.Kx,dom.Ky,gradeta2Lphi_hat);
            Gradphi_hat    =funOprt_grad2d_hat(dom.Kx,dom.Ky,phi_hat);
            DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,Gradphi_hat);
            DivGradphi     =funC_ifft2(DivGradphi_hat);
            Leta2DivGradphi_hat=Oprt.L2d.*fft2(eta.^2.*DivGradphi);
            delphiH3_hat  =LetaLetaLphi_hat+0.5.*(DivGradeta2Lphi_hat+Leta2DivGradphi_hat);
            deletaH3_hat  =fft2(Lphi.*(eta.*DivGradphi+LetaLphi));
        end
    end
    delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat+delphiH3_hat).*ChiAdj);
    deletaH_hat=deletaH1_hat+fft2(funC_ifft2(deletaH2_hat+deletaH3_hat).*ChiAdj);
end

if strcmp(model.breaking.check,'Yes')
     [B] = breaking_process2d(eta,gradphi,dom,parBreak,time,g,t_init);
     Sb_hat  = -parBreak.delb.^2*fft2(B.*funC_ifft2(delphiH_hat).^2);%parBreak.delb
   
%     H=(Depth+eta);
%     [B] = breaking_process2d(eta,gradphi,dom,parBreak,time,g);
%     visc=B.*parBreak.delb.^2.*(H).*funC_ifft2(delphiH_hat);
%     Flux=-visc.*funC_ifft2(delphiH_hat);Kx=dom.Kx;
%     [gradFlux_x,gradFlux_y]=gradient(Flux,dom.dx,dom.dy);
%     Sb_hat=fft2((cumtrapz(dom.X',gradFlux_x'./H'))');   
    
else
    Sb_hat=0;
end

if dom.friction.input.check==1
    if model.nonlinear==1
        gradphi       = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
    end
    absGradphi=sqrt(gradphi.x.^2+gradphi.y.^2);
    Depth  =-dom.bathy.profile ;
    Sf=-dom.friction.Char.*(1./(Depth+eta).^(4/3)).*phi.*absGradphi;
    Sf_hat=fft2(Sf);
else
    Sf_hat=0;
end

dteta_hat=(delphiH_hat+Source_hat-dampcoef*fft2(eta.*dampchar)).*aal;
dtphi_hat=(-deletaH_hat+Sb_hat+Sf_hat-dampcoef*fft2(phi.*dampchar)).*aal;

dz1_hat=reshape(dteta_hat,1,[]);
dz2_hat=reshape(dtphi_hat,1,[]);

dz=[dz1_hat,dz2_hat].';


if par.IP.check==1
    if time>=par.IP.time(1) && time<=par.IP.time(2)
        if time>=par.IP.time(1)+(iterInteriordt-1)*dt*par.IP.time(3)
            dteta(iterInterior,:,:)   =funC_ifft2(dteta_hat);
            dtphihat(iterInterior,:,:)=dtphi_hat;
            timeIP(iterInterior,1)    =time;
            phiIP(iterInterior,:,:)   =phi;
            etaIP(iterInterior,:,:)   =eta;
            
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end



if plotchecker==1 
[X,Y]=meshgrid(dom.X,dom.Y);
if  iterdtcheck==1
    f_check=figure;
    set(f_check,'Renderer','zbuffer'); 
end

    if mod(iterdtcheck,2)
       surf(X,Y,Refl_Eta,'Edgecolor','none')
        view(2)
        colorbar;
        pause(0.001);
    end

 iterdtcheck= iterdtcheck+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if max(any(abs(eta)>6*par.Hs))==1 %%
%     Idstop=2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
