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

function [dz] = RHSFlat2dWalluv(time,z_hat,input,model,par,parBreak,dom,influx,ivp,...
    Oprt,timeSimul,ProgStatbar)

plotchecker=0;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=dom.fbl.char; %dampchar=char./length of fourier boundary;
dampcoef=7*Oprt.Cpeak;
aal     =dom.aal;
Depth  =-dom.bathy.profile ;
ChiAdj =influx.ChiAdj;
Wallchar=dom.wall.char;
GradWallchar=dom.wall.gradchar;
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
z1_hat=z_hat(1:N);     eta_hat=reshape(z1_hat,dom.Ny,[]);
z2_hat=z_hat(N+1:2*N); u_hat=reshape(z2_hat,dom.Ny,[]);
z3_hat=z_hat(2*N+1:end); v_hat=reshape(z3_hat,dom.Ny,[]);
eta=funC_ifft2(eta_hat);u=funC_ifft2(u_hat);v=funC_ifft2(v_hat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gradphi.x          = u;gradphi.y       = v;
gradphi_hat.x_hat  = u_hat;gradphi_hat.y_hat = v_hat;
H=Depth+eta;
%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model.evol,'SWE')
    
    if model.nonlinear==1
        depthgradphi_hat.x_hat=fft2(funC_ifft2(Depth.*gradphi_hat.x_hat).*Wallchar);
        depthgradphi_hat.y_hat=fft2(funC_ifft2(Depth.*gradphi_hat.y_hat).*Wallchar);
        div_depthgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,depthgradphi_hat);
        delphiH_hat  = -div_depthgradphi_hat;
        deletaH_hat  = g.*eta_hat;
    else
        Hgradphi_hat.x_hat=fft2(H.*gradphi.x.*Wallchar);
        Hgradphi_hat.y_hat=fft2(H.*gradphi.y.*Wallchar);
        div_Hgradphi_hat= funOprt_div2d(dom.Kx,dom.Ky,Hgradphi_hat);
        gradphi2      = funOprt_innerproduct(gradphi,gradphi);
        delphiH_hat  = -div_Hgradphi_hat;
        deletaH_hat  = g.*eta_hat+0.5.*fft2(gradphi2.*Wallchar);
    end
    dx_deletaH_hat =1i.*dom.Kx.*deletaH_hat;
    dy_deletaH_hat =1i.*dom.Ky.*deletaH_hat;
    
elseif strcmp(model.evol,'HS')
    if model.nonlinear==1
        if strcmp(dom.wall.ReflCoef.AllOneFlagCutK,'True')
            Cp2gradphi_hat.x_hat = Oprt.Cp2d.^2.*fft2(gradphi.x);
            Cp2gradphi_hat.y_hat = Oprt.Cp2d.^2.*fft2(gradphi.y);
            
            Flux.x               = funC_ifft2(Cp2gradphi_hat.x_hat)/g;
            Flux.y               = funC_ifft2(Cp2gradphi_hat.y_hat)/g;
            M0uvwall_hat= fft2(funC_ifft2(-funOprt_div2d(dom.Kx,dom.Ky,Cp2gradphi_hat)/g).*Wallchar)...
                - fft2(funOprt_innerproduct(Flux,GradWallchar));
        else
            Cu.x                 =funC_ifft2(Oprt.Cp2d.*gradphi_hat.x_hat);
            Cu.y                 =funC_ifft2(Oprt.Cp2d.*gradphi_hat.y_hat);
            Cp2gradphiWall_hat.x_hat = Oprt.Cp2d.*fft2(Cu.x.*Wallchar);
            Cp2gradphiWall_hat.y_hat = Oprt.Cp2d.*fft2(Cu.y.*Wallchar);
            
            M0uvwall_hat= -funOprt_div2d(dom.Kx,dom.Ky,Cp2gradphiWall_hat)/g;
        end
        %%%%%% linear term %%%%%%%%%%%
        delphiH1_hat  =M0uvwall_hat;
        deletaH1_hat  =g.*eta_hat;
        delphiH2_hat  =0; deletaH2_hat = 0;
        delphiH3_hat  =0; deletaH3_hat=0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        Cu.x                 =funC_ifft2(Oprt.Cp2d.*gradphi_hat.x_hat);
        Cu.y                 =funC_ifft2(Oprt.Cp2d.*gradphi_hat.y_hat);
        Cp2gradphiWall_hat.x_hat = Oprt.Cp2d.*fft2(Cu.x.*Wallchar);
        Cp2gradphiWall_hat.y_hat = Oprt.Cp2d.*fft2(Cu.y.*Wallchar);
        
        M0uvwall_hat= -funOprt_div2d(dom.Kx,dom.Ky,Cp2gradphiWall_hat)/g;
        %%%%%% linear term %%%%%%%%%%%
        delphiH1_hat  =M0uvwall_hat;
        deletaH1_hat  =g.*eta_hat;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         C2u_hat.x_hat=Oprt.Cp2dsq.^2.*gradphi_hat.x_hat;
%         C2u_hat.y_hat=Oprt.Cp2dsq.^2.*gradphi_hat.y_hat;
        
        M0uv              = funC_ifft2(M0uvwall_hat);%-funC_ifft2(funOprt_div2d(dom.Kx,dom.Ky,C2u_hat))/g;
        LetaM0uvWall_hat  = Oprt.L2d.*fft2(eta.*M0uv.*Wallchar);
        etagradphiWall_hat.x_hat=fft2(eta.*gradphi.x.*Wallchar);
        etagradphiWall_hat.y_hat=fft2(eta.*gradphi.y.*Wallchar);
        divetagradphiWall_hat = funOprt_div2d(dom.Kx,dom.Ky,etagradphiWall_hat);
        gradphi2  = funOprt_innerproduct(gradphi,gradphi);
        delphiH2_hat  = -(LetaM0uvWall_hat+divetagradphiWall_hat);
        deletaH2_hat  = 0.5*fft2((gradphi2-M0uv.^2).*Wallchar);
        delphiH3_hat  =0; deletaH3_hat=0;
        
        if model.nonlinear~=2
            LetaM0uv        =funC_ifft2(Oprt.L2d.*fft2(eta.*M0uv));
            LetaWallLetaLphi_hat=Oprt.L2d.*fft2(eta.*Wallchar.*LetaM0uv);
            LetaLetaWallLphi_hat=Oprt.L2d.*fft2(eta.*funC_ifft2(LetaM0uvWall_hat));
            
            Walleta2M0uv_hat=fft2(Wallchar.*eta.^2.*M0uv);
            gradWalleta2M0uv_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,Walleta2M0uv_hat);
            DivGradWalleta2M0uv_hat=funOprt_div2d(dom.Kx,dom.Ky,gradWalleta2M0uv_hat);
            DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,gradphi_hat);
            DivGradphi     =funC_ifft2(DivGradphi_hat);
            Leta2WallDivGradphi_hat=Oprt.L2d.*fft2(eta.^2.*Wallchar.*DivGradphi);
            delphiH3_hat  =0.5*(LetaWallLetaLphi_hat+LetaLetaWallLphi_hat+...
                Leta2WallDivGradphi_hat+DivGradWalleta2M0uv_hat);
            deletaH3_hat  =fft2(M0uv.*(eta.*DivGradphi.*Wallchar+...
                0.5*funC_ifft2(LetaM0uvWall_hat)+0.5*Wallchar.*LetaM0uv));
        end
    end
    delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat+delphiH3_hat).*ChiAdj);
    dx_deletaH_hat=1i.*dom.Kx.*deletaH1_hat+fft2(funC_ifft2(1i.*dom.Kx.*(deletaH2_hat+deletaH3_hat)).*ChiAdj);
    dy_deletaH_hat=1i.*dom.Ky.*deletaH1_hat+fft2(funC_ifft2(1i.*dom.Ky.*(deletaH2_hat+deletaH3_hat)).*ChiAdj);
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

if dom.wall.NInfl>0
    [Swall_hat]=funW_influxsource_2D(influx,dom,eta,time);
else
    Swall_hat=0;
end

dteta_hat=(delphiH_hat+Source_hat+Swall_hat-dampcoef*fft2(eta.*dampchar)).*aal;
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
    
    if mod( iterdtcheck,2)
        % subplot(2,1,1)
        surf(X,Y,funC_ifft2(M0uvwall_hat),'Edgecolor','none')
        view(2)
        colorbar;
        %    caxis([-0.1;0.1]);
        % zlim([-0.1;0.1])
        title(['time= ',num2str(time)])
        %         subplot(2,1,2)
        %         plot(dom.Y,INsig_t0,'r',dom.Y,eta(:,funC_closest(dom.X,60)),'--b',...
        %             dom.Y,Source(:,funC_closest(dom.X,60)),'-.g')
        plot_properties;
        % ylim([-0.1 0.1])
        %caxis([-0.1;0.1])
        pause(0.001);
    end
    iterdtcheck= iterdtcheck+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if max(any(abs(eta)>6*par.Hs))==1 %%
%     Idstop=2;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

