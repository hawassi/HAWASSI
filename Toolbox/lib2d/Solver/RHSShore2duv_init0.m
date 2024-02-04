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
%%%%%%%%%% version: 24 March 2020                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% RHS2D for shore/ runup                                 %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dz] = RHSShore2duv_init0(time,z_hat,input,model,par,bath,parBreak,dom,influx,ivp,...
    Oprt,timeSimul,ProgStatbar)
plotchecker=0;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=1-dom.cfSA;
dampcoef=Oprt.Cpeak;
aal     =dom.aal;
Depth   =-dom.bathy.profile ;
ChiAdj  =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
H_minShore  = dom.bathy.HminShore;
H_minDisp   = dom.bathy.HminDisp;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar iterdtcheck
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
eta=funC_ifft2(eta_hat);eta_prev=eta;
u=funC_ifft2(u_hat);
v=funC_ifft2(v_hat);

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;

%%%Define active domain by a characteristic function
WaveChar=ones(size(eta));

H=ChiAdj.*eta-dom.bathy.min;
WaveChar(H<H_minShore)=zeros;

if strcmp(bath.type,'Shore (in x-axis)') || strcmp(bath.type,'Userdefined (shore in x-axis)')
%     IndxZerofirst=funC_findIndexinColumn(WaveChar,'first');
%     IndxZerolast =funC_findIndexinColumn(WaveChar,'last');
    for jj=1:dom.Ny
        IndxZerofirst=find(WaveChar(jj,:) == 0, 1,'first');
        IndxZerolast=find(WaveChar(jj,:) == 0, 1,'last');
        WaveChar(jj,IndxZerofirst:IndxZerolast)=0;
         %%% added smoother damp in the shore, 20201203
         dampchar(jj,IndxZerofirst:IndxZerolast)=1;
         if IndxZerofirst==1 %% propagation to the left
             dampchar(jj,IndxZerofirst:IndxZerolast+1)=smooth(dampchar(jj,IndxZerofirst:IndxZerolast+1),3);
         else
             dampchar(jj,IndxZerofirst-1:IndxZerolast)=smooth(dampchar(jj,IndxZerofirst-1:IndxZerolast),3);
         end
    end
   

else
%     IndxZerofirst=funC_findIndexinColumn(WaveChar.','first');
%     IndxZerolast =funC_findIndexinColumn(WaveChar.','last');
    for jj=1:dom.Nx
        IndxZerofirst=find(WaveChar(:,jj) == 0, 1,'first');
        IndxZerolast=find(WaveChar(:,jj) == 0, 1,'last');
        WaveChar(IndxZerofirst:IndxZerolast,jj)=0;
        %%% added smoother damp in the shore, 20201203
        dampchar(IndxZerofirst:IndxZerolast,jj)=1;
        if IndxZerofirst==1 %% propagation to the bottom
            dampchar(IndxZerofirst:IndxZerolast+1,jj)=smooth(dampchar(IndxZerofirst:IndxZerolast+1,jj),3);
        else
            dampchar(IndxZerofirst-1:IndxZerolast,jj)=smooth(dampchar(IndxZerofirst-1:IndxZerolast,jj),3);
        end
    end
end

Wallchar=dom.wall.char;
ShoreChar=1;
HeavChar=ShoreChar.*Wallchar;
gradphi.x = u;gradphi.y = v;
gradphi_hat.x_hat = u_hat;gradphi_hat.y_hat = v_hat;

absuv2    = funOprt_innerproduct(gradphi,gradphi);

if strcmpi(model.runupevo,'Hsdirect')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Cp2dsqInterp = Oprt.Cp2dSqInterp;
    CharH  = ones(size(eta));
    CharH(H<H_minDisp)= 0;
    Cp2dsqInterp.gammin1  = ppval(Cp2dsqInterp.gammin1,H).*CharH;
    Cp2dsqInterp.gamplus1 = ppval(Cp2dsqInterp.gamplus1,H).*CharH;
    Cp2dsqInterp.gammid1  = ppval(Cp2dsqInterp.gammid1,H).*CharH;
    Cp2dsqInterp.gammin2  = ppval(Cp2dsqInterp.gammin2,H).*CharH;
    Cp2dsqInterp.gamplus2 = ppval(Cp2dsqInterp.gamplus2,H).*CharH;
    Cp2dsqInterp.gammid2  = ppval(Cp2dsqInterp.gammid2,H).*CharH;
    
    C2gradphi_hat.x_hat = funOprt_Cp2d_runup(Cp2dsqInterp,u_hat,u);
    C2gradphi_hat.y_hat = funOprt_Cp2d_runup(Cp2dsqInterp,v_hat,v);
    
    C2gradphi.x = funC_ifft2(C2gradphi_hat.x_hat);
    C2gradphi.y = funC_ifft2(C2gradphi_hat.y_hat);
    
    C2gradphi.x(H<H_minDisp) = g*H(H<H_minDisp).*u(H<H_minDisp);
    C2gradphi.y(H<H_minDisp) = g*H(H<H_minDisp).*v(H<H_minDisp);
    
    C2gradphi.x = C2gradphi.x.*Wallchar;
    C2gradphi.y = C2gradphi.y.*Wallchar;
    
    C2gradphi_hat.x_hat = fft2(C2gradphi.x);
    C2gradphi_hat.y_hat = fft2(C2gradphi.y);
    
    Lphi_hat = -1/g*funOprt_div2d(dom.Kx,dom.Ky,C2gradphi_hat);
    Lphi = funC_ifft2(Lphi_hat);
    
    gEta=g.*(eta+dom.bathy.plus);
    deletaH_hat = fft2(gEta) + 1/2*fft2((absuv2-Lphi.^2).*ChiAdj.*Wallchar);
    delphiH_hat =-funOprt_div2d(dom.Kx,dom.Ky,C2gradphi_hat)/g;
    
    dx_deletaH_hat=1i.*dom.Kx.*deletaH_hat;
    dy_deletaH_hat=1i.*dom.Ky.*deletaH_hat;
else
    
    Cp2dsqInterp = Oprt.Cp2dSqInterp; Cp2dInterp=Oprt.Cp2dInterp;Om2dsqInterp=Oprt.Om2dSqInterp;
    
    Cpgradphi.x   =funC_ifft2(funOprt_Cp2d_runup(Cp2dInterp, gradphi_hat.x_hat,gradphi.x));
    Cpgradphi.y   =funC_ifft2(funOprt_Cp2d_runup(Cp2dInterp, gradphi_hat.y_hat,gradphi.y));
    CgradphiChar.x     =Cpgradphi.x.*HeavChar;
    CgradphiChar.y     =Cpgradphi.y.*HeavChar;
    
    CgradphiChar_hat.x_hat = fft2(CgradphiChar.x);
    CgradphiChar_hat.y_hat = fft2(CgradphiChar.y);
    C2gradphiChar_hat.x_hat = funOprt_Cp2d_runup(Cp2dInterp,CgradphiChar_hat.x_hat,CgradphiChar.x);
    C2gradphiChar_hat.y_hat = funOprt_Cp2d_runup(Cp2dInterp,CgradphiChar_hat.y_hat,CgradphiChar.y);
    
    delphiH1_hat=-1/g*funOprt_div2d(dom.Kx,dom.Ky,C2gradphiChar_hat);
    gEta_hat=fft2(g.*(eta+dom.bathy.plus));
    deletaH1_hat  =gEta_hat;
    
    eta(Depth<H_minDisp)   = H(Depth<H_minDisp);
    
    C2gradphi_hat.x_hat    = funOprt_Cp2d_runup(Cp2dsqInterp,u_hat,u);
    C2gradphi_hat.y_hat    = funOprt_Cp2d_runup(Cp2dsqInterp,v_hat,v);
    M0uv_hat               = -1/g*funOprt_div2d(dom.Kx,dom.Ky,C2gradphi_hat);
    M0uv                   = funC_ifft2(M0uv_hat);
    etaM0uvChar             = eta.*M0uv.*HeavChar;
    etaM0uvChar_hat         = fft2(etaM0uvChar);
    LetaM0uvChar_hat        = funOprt_L2d_runup(g,Om2dsqInterp,etaM0uvChar_hat,etaM0uvChar);
    etagradphi_hat.x_hat   = fft2(eta.*gradphi.x);
    etagradphi_hat.y_hat   = fft2(eta.*gradphi.y);
    divetagradphi_hat      = funOprt_div2d(dom.Kx,dom.Ky,etagradphi_hat);
    delphiH2_hat           = -(LetaM0uvChar_hat+divetagradphi_hat);
    deletaH2_hat           = 0.5.*fft2((absuv2-M0uv.^2).*HeavChar);
    if model.nonlinear==2
        delphiH3_hat           = 0;deletaH3_hat=0;
    else
        etaM0uv             =eta.*M0uv;
        LetaM0uv_hat        =funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaM0uv),etaM0uv);
        LetaM0uv            =funC_ifft2(LetaM0uv_hat);
        etaCharLetaM0uv     =eta.*HeavChar.*LetaM0uv;
        LetaCharLetaM0uv_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaCharLetaM0uv),etaCharLetaM0uv);
        etaLetaM0uvChar     =eta.*funC_ifft2(LetaM0uvChar_hat);
        LetaLetaCharM0uv_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(etaLetaM0uvChar),etaLetaM0uvChar);
        Chareta2M0uv_hat=fft2(Wallchar.*eta.^2.*M0uv);
        gradChareta2M0uv_hat=funOprt_grad2d_hat(dom.Kx,dom.Ky,Chareta2M0uv_hat);
        DivGradChareta2M0uv_hat=funOprt_div2d(dom.Kx,dom.Ky,gradChareta2M0uv_hat);
        DivGradphi_hat =funOprt_div2d(dom.Kx,dom.Ky,gradphi_hat);
        DivGradphi     =funC_ifft2(DivGradphi_hat);
        eta2CharDivGradphi=eta.^2.*HeavChar.*DivGradphi;
        Leta2CharDivGradphi_hat=funOprt_L2d_runup(g,Om2dsqInterp,fft2(eta2CharDivGradphi),eta2CharDivGradphi);
        delphiH3_hat  =0.5*(LetaCharLetaM0uv_hat+LetaLetaCharM0uv_hat+...
            Leta2CharDivGradphi_hat+DivGradChareta2M0uv_hat);
        deletaH3_hat  =fft2(M0uv.*(eta.*DivGradphi.*HeavChar+...
            0.5*funC_ifft2(LetaM0uvChar_hat)+0.5*HeavChar.*LetaM0uv));
    end
    delphiH_hat=delphiH1_hat+fft2(funC_ifft2(delphiH2_hat+delphiH3_hat).*ChiAdj);
    dx_deletaH_hat=1i.*dom.Kx.*deletaH1_hat+fft2(funC_ifft2(1i.*dom.Kx.*(deletaH2_hat+deletaH3_hat)).*ChiAdj);
    dy_deletaH_hat=1i.*dom.Ky.*deletaH1_hat+fft2(funC_ifft2(1i.*dom.Ky.*(deletaH2_hat+deletaH3_hat)).*ChiAdj);
end

if strcmp(model.breaking.check,'Yes')
    parBreak.char=ones(size(dom.XX));
    parBreak.char(WaveChar==0)=0;%update
    [B] = breaking_process2d(eta_prev,gradphi,dom,parBreak,time,g,t_init);
    visc=B.*parBreak.delb.^2.*H.*funC_ifft2(delphiH_hat);
    Flux=-visc.*funC_ifft2(delphiH_hat);
    [dxFlux, dyFlux]=gradient(Flux,dom.dx,dom.dy);
    oneperH=1./H;oneperH(WaveChar==0)=0;
    Sb_hat.x  = fft2(dxFlux.*oneperH);
    Sb_hat.y  = fft2(dyFlux.*oneperH);
else
    Sb_hat.x=0; Sb_hat.y=0;
end

if dom.friction.input.check==1
    absGradphi=sqrt(gradphi.x.^2+gradphi.y.^2);
    oneperH43=1./(H.^(4/3));oneperH43(WaveChar==0)=0;
    Sfx=-dom.friction.Char.*oneperH43.*u.*absGradphi;
    Sfy=-dom.friction.Char.*oneperH43.*v.*absGradphi;
    Sf_hat.x=fft2(Sfx);
    Sf_hat.y=fft2(Sfy);
else, Sf_hat.x=0;Sf_hat.y=0;
end


Upeta=dampcoef*eta_prev; Upu=dampcoef*u; Upv=dampcoef*v;
Upeta(WaveChar==0)=dampcoef*H(WaveChar==0);

dteta_hat=(delphiH_hat+Source_hat-fft2(Upeta.*dampchar)).*aal;
dtu_hat=(-dx_deletaH_hat+Sb_hat.x+Sf_hat.x-fft2(Upu.*dampchar)).*aal;
dtv_hat=(-dy_deletaH_hat+Sb_hat.y+Sf_hat.y-fft2(Upv.*dampchar)).*aal;

% figure;
% plot(eta(1,:));
% figure
% plot(phi(1,:));
% pause
%
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
            etaIP(iterInterior,:,:)   =eta_prev+dom.bathy.plus;
            iterInterior=iterInterior+1;
            iterInteriordt=iterInteriordt+1;
        end
    end
end

if plotchecker==1 && time >800
    stepS=1;
    [X,Y]=meshgrid(dom.X(1:stepS:end),dom.Y(1:stepS:end));
    if  iterdtcheck==1
        figure;
    end
    if mod( iterdtcheck,2)
    
        surf(X,Y,eta(1:stepS:end,1:stepS:end),'Edgecolor','none')
        view(2)
        colorbar;
        axis([1000 2200 0 100]);
        caxis([-1;1]);
%         zlim([-0.05;0.05])
        title(['time= ',num2str(time)])
%         pause(0.01);
    end
    iterdtcheck= iterdtcheck+1;
end
