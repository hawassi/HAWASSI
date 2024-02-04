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
function [dz] = RHSShore2d_init0(time,z_hat,input,model,par,bath,parBreak,dom,influx,ivp,...
                Oprt,timeSimul,ProgStatbar)
plotchecker=1;
%%%%%passing parameter%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=par.g;
dampchar=1-dom.cfSA; 
dampcoef=Oprt.Cpeak;
aal     =dom.aal;
Depth   =-dom.bathy.profile ;
ChiAdj  =influx.ChiAdj;
jProgressBar=ProgStatbar.jProgressBar;
statusbarObj=ProgStatbar.statusbarObj;
H_minShore= dom.bathy.HminShore;
H_min     = 5*H_minShore;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global flagWarnIVP iterProgBar 
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
eta=funC_ifft2(eta_hat);
phi=funC_ifft2(phi_hat);

%%update phi
% phi = phi - mean(mean(phi));
% phi_hat= fft2(phi);

%%%%%%%%%%%%%%%%%%Influxing%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Wave_Generation_ODE;

%%%Define active domain by a characteristic function
WaveChar=ones(size(eta));
H=ChiAdj.*eta+Depth;
H(Depth<0)=ChiAdj(Depth<0).*eta(Depth<0);
WaveChar(H<H_minShore)=zeros;

if strcmp(bath.type,'Shore (in x-axis)') %corrected nida
    IndxZerofirst=funC_findIndexinColumn(WaveChar,'first');
    IndxZerolast =funC_findIndexinColumn(WaveChar,'last');
    WaveChar(:,IndxZerofirst(:):IndxZerolast(:))=0;   
else
    IndxZerofirst=funC_findIndexinColumn(WaveChar','first');
    IndxZerolast =funC_findIndexinColumn(WaveChar','last');
    WaveChar(IndxZerofirst(:):IndxZerolast(:),:)=0;
end

dampchar(WaveChar==0)=ones;
Wallchar=dom.wall.char;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Cp2dsq = Oprt.Cp2dSqInterp;
CharH  = ones(size(eta));
CharH(H<H_min)= 0;
Cp2dsq.gammin1  = ppval(Cp2dsq.gammin1,H).*CharH; 
Cp2dsq.gamplus1 = ppval(Cp2dsq.gamplus1,H).*CharH; 
Cp2dsq.gammid1  = ppval(Cp2dsq.gammid1,H).*CharH; 
Cp2dsq.gammin2  = ppval(Cp2dsq.gammin2,H).*CharH; 
Cp2dsq.gamplus2 = ppval(Cp2dsq.gamplus2,H).*CharH; 
Cp2dsq.gammid2  = ppval(Cp2dsq.gammid2,H).*CharH; 

%%phi to u
gradphi = funOprt_grad2d(dom.Kx,dom.Ky,phi_hat);
gradphi_hat.x_hat = fft2(gradphi.x);
gradphi_hat.y_hat = fft2(gradphi.y);

%% general
C2gradphi_hat.x_hat = funOprt_Cp2d_runup(Cp2dsq,gradphi_hat.x_hat,gradphi.x); 
C2gradphi_hat.y_hat = funOprt_Cp2d_runup(Cp2dsq,gradphi_hat.y_hat,gradphi.y); 

C2gradphi.x = funC_ifft2(C2gradphi_hat.x_hat);
C2gradphi.y = funC_ifft2(C2gradphi_hat.y_hat);

C2gradphi.x(H<H_min) = g*H(H<H_min).*gradphi.x(H<H_min); 
C2gradphi.y(H<H_min) = g*H(H<H_min).*gradphi.y(H<H_min); 

C2gradphi.x = C2gradphi.x.*Wallchar;
C2gradphi.y = C2gradphi.y.*Wallchar;

C2gradphi_hat.x_hat = fft2(C2gradphi.x); 
C2gradphi_hat.y_hat = fft2(C2gradphi.y); 

Lphi_hat = -1/g*funOprt_div2d(dom.Kx,dom.Ky,C2gradphi_hat);
Lphi = funC_ifft2(Lphi_hat);

gradphi2 = funOprt_innerproduct(gradphi,gradphi);

gEta=g.*eta;
gEta(Depth<0)=g.*(eta(Depth<0)-Depth(Depth<0));

delphiH_hat = Lphi_hat;
deletaH_hat = fft2(gEta) + 1/2*fft2((gradphi2-Lphi.^2).*ChiAdj.*Wallchar); 

%% Longcrested case in Shore-x
% CCu_hat = funOprt_Cp2d_runup(Cp2dsq,gradphi_hat.x_hat,gradphi.x);  
% CCu_ShoreChar  = funC_ifft2(CCu_hat);  
%   
% CCu_ShoreChar(H<H_min) = g*H(H<H_min).*gradphi.x(H<H_min);    
% CCu_ShoreChar_hat      = fft2(CCu_ShoreChar);    
%  
% kC2_u = funC_ifft2(dom.Kx.*CCu_ShoreChar_hat);
% 
% delphiH_hat = -1i.*dom.Kx.*1/g.*CCu_ShoreChar_hat;
% deletaH_hat = (fft2(g.*eta) + fft2(1/2*funC_ifft2(fft2(u.^2-(kC2_u./g).^2)).*ChiAdj.*Wallchar)); 

if strcmp(model.breaking.check,'Yes')
    parBreak.char(WaveChar==0)=0;%update 
    [B] = breaking_process2d(eta,gradphi,dom,parBreak,time,g);
    Sb_hat  = -parBreak.delb.^2*fft2(B.*funC_ifft2(delphiH_hat).^2);
else
    Sb_hat=0;
end

if dom.friction.input.check==1
    Sf=-dom.friction.Char.*(1./H).*phi;
    Sf(H<H_min)=0;
    Sf_hat=fft2(Sf);
else
    Sf_hat=0;
end

damptype = 1;
if (damptype == 0)%original
    %%Potential formulation
    
    Upeta=dampcoef*eta;
    Upphi=dampcoef*phi;
    Upeta(WaveChar==0)=dampcoef*H(WaveChar==0);%%update dry area
    dteta_hat=(delphiH_hat+Source_hat-fft2(Upeta.*dampchar)).*aal;
    dtphi_hat=(-deletaH_hat+Sb_hat+Sf_hat-fft2(Upphi.*dampchar)).*aal;
else
    %%Velocity formulation
    if strcmp(bath.type,'Shore (in x-axis)')
        Upeta=dampcoef*eta;
        Upeta(WaveChar==0)=dampcoef*H(WaveChar==0);%%update dry area

        Updxphi=dampcoef*gradphi.x;
        Updyphi=dampcoef*gradphi.y;
        Upphi=dampcoef*phi;
        graddampchar = funOprt_grad2d(dom.Kx,dom.Ky,fft2(dampchar));
       
        dteta_hat=(delphiH_hat+Source_hat-fft2(Upeta.*dampchar)).*aal;
 %       corrdamp1_hat=fft2(graddampchar.x.*Upphi)./(1i*dom.Kx);
%        corrdamp1=funC_ifft2(corrdamp1_hat);
%        corrdamp1=corrdamp1-mean(mean(corrdamp1));
%        corrdamp1_hat=fft(corrdamp1);
%        corrdamp1_hat(1,:)=0;
%        corrdamp1_hat(:,1)=0;
       
 %          corrdamp2_hat=(fft2(graddampchar.y.*Updxphi-graddampchar.x.*Updyphi)./(1i*dom.Kx))./(1i*dom.Ky); 
         corrdamp2_hat= fft2(funC_ifft2(fft2(graddampchar.y.*Updxphi)./(1i*dom.Kx))...
                            -funC_ifft2(fft2(graddampchar.x.*Updyphi)./(1i*dom.Kx)))./(1i*dom.Ky); 
        corrdamp2_hat(1,:)=0;
        corrdamp2_hat(:,1)=0;

        damp0_hat=fft2(Updxphi.*dampchar)./(1i*dom.Kx);
        % damp0_hat=-fft2(cumtrapz(dom.X,Updxphi.*dampchar,2));
        % damp1_hat=-fft2(Upphi.*dampchar)+corrdamp1_hat;
        dtphi_hat=(-deletaH_hat+Sb_hat+Sf_hat-(damp0_hat- corrdamp2_hat)).*aal;
     
    else
        Upeta=dampcoef*eta;
        Upphi=dampcoef*gradphi.y;
        Upeta(WaveChar==0)=dampcoef*H(WaveChar==0);%%update dry area
        Upphi(WaveChar==0)=dampcoef*gradphi.y(WaveChar==0);%%update dry area
        Cphi = ones(size(phi))*mean(mean(phi));
        Cphihat = 0;%fft2(Cphi);

        dteta_hat=(delphiH_hat+Source_hat-fft2(Upeta.*dampchar)).*aal;
        dtphi_hat=(-deletaH_hat+Sb_hat+Sf_hat-fft2(Upphi.*dampchar)./(1i*dom.Ky)-Cphihat).*aal;
    end
end

% figure;
% plot(eta(1,:));
% figure
% plot(phi(1,:));
% pause

dz1_hat=reshape(dteta_hat,1,[]);
dz2_hat=reshape(dtphi_hat,1,[]);
dz=[dz1_hat,dz2_hat].';

% if plotchecker==1
%     stepS=4;
%     [X,Y]=meshgrid(dom.X(1:stepS:end),dom.Y(1:stepS:end));
%     if  iterdtcheck==1
%         figure;
%     end
%     if mod( iterdtcheck,2)
%         hold off;
%         plot(dom.X(1,:),eta(1,:))         
% %         surf(X,Y,eta(1:stepS:end,1:stepS:end),'Edgecolor','none')
%         % hold on
%         % surf(X,Y,WaveChar(1:stepS:end,1:stepS:end),'Edgecolor','none')
%         view(2)
%         colorbar;
%         axis([-3500 100 -0.3 0.3]);
%         caxis([-0.1;0.1]);
% %         zlim([-0.05;0.05])
%         title(['time= ',num2str(time)])
% %         pause(0.01);
%     end
%     iterdtcheck= iterdtcheck+1;
% end
