%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% HAWASSI-AB 2D                                         %%%%%%%%%%
%%%%%%%%%% Hamiltonian Wave-Ship-Structure Interaction           %%%%%%%%%%
%%%%%%%%%% Copyright (c): LabMath-Indonesia                      %%%%%%%%%%
%%%%%%%%%% version: 5 July 2016                                  %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%RK%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Short-crested wave generator                          %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [INsig,SpectProp]=funG_JONSWAP_2D_signal(g,dom,inputWave,timesig,halfomsig,spatinflux,depthI,I)
gam=inputWave.JS_gamma(I);
Tp=inputWave.Tp(I);
Hs=inputWave.Hs(I);
JS_s=inputWave.JS_s(I);
mDir_rad=inputWave.direction(I)*pi./180;

domg=halfomsig(2)-halfomsig(1);
Nomg=length(halfomsig);

%%Preparation of 1D-JONSWAP-spectrum
[JS_spect1d]=funG_JONSWAP_1D_spectrum(g,gam,halfomsig,Tp,Hs);

%%Preparation of directional spreading
if JS_s==inf  % for LONG-CRESTED WAVE
    Rtheta = zeros(1,Nomg);                  % uniform for long-crested waves
    theta  = linspace(-pi/2,pi/2, Nomg ); 
    StdTheta = 0;cos2pdf=zeros(size(halfomsig));
    cos2pdf(funC_closest(theta,0))=1/(theta(2)-theta(1));
else                  % for SHORT-CRESTED WAVES
    [Rtheta, StdTheta,cos2pdf,theta] = funG_dir_spreading_cos2s(Nomg, JS_s,mDir_rad);
end

SpectProp.varDensSpect=JS_spect1d;
SpectProp.StdTheta=StdTheta;
SpectProp.cos2pdf=cos2pdf;
SpectProp.mDir=mDir_rad;
SpectProp.halfomsig=halfomsig;
a_bar    = sqrt(2*JS_spect1d*domg);
rng('shuffle');pause(0.1); %%shuffle random generator; 
RP = rand(Nomg,1);                  % random phase : the same size as omg and theta
RP = RP.*(2*pi);                    % phase \in [0,2pi]

XYinfl=spatinflux.line(I).xy;

if strcmp(spatinflux.line(I).Orientation,'Vertical')
 INsig=zeros(Nomg,dom.Ny);
 Indymin=funC_closest(dom.Y,min(XYinfl(:,2)));
 Indymax=funC_closest(dom.Y,max(XYinfl(:,2)));
 
 Depth_influx_line=zeros(dom.Ny,1);
 Depth_influx_line(Indymin:Indymax)=depthI;
 Depth_influx_line(1:Indymin-1)=depthI(1);
 Depth_influx_line(Indymax+1:end)=depthI(end);
 Xline=zeros(dom.Ny,1);   
 Xline(Indymin:Indymax)=XYinfl(:,1);
 Xline(1:Indymin-1)=XYinfl(1,1);
 Xline(Indymax+1:end)=XYinfl(end,1);
 Xline=repmat(Xline',Nomg,[]);
 
 [Y,WW]=meshgrid(dom.Y,halfomsig);
 KK=funOprt_invOmExact(dom.ky,halfomsig,Depth_influx_line,g);
 Rtheta=repmat(Rtheta',1,dom.Ny);
 RP=repmat(RP,1,dom.Ny);
 a_bar=repmat(a_bar',1,dom.Ny);
 
  for ii=1:length(timesig)          
      eta_y_kjj_thetakk=a_bar.*cos(KK.*sin(Rtheta+mDir_rad).*Y...
                                +KK.*cos(Rtheta+mDir_rad).*Xline...
                                 -WW.*timesig(ii)+RP);
      INsig(ii,:)=sum(eta_y_kjj_thetakk);
  end
else
 INsig=zeros(Nomg,dom.Nx);
 Indxmin=funC_closest(dom.X,min(XYinfl(:,1)));
 Indxmax=funC_closest(dom.X,max(XYinfl(:,1)));
 
 Depth_influx_line=zeros(dom.Nx,1);
 Depth_influx_line(Indxmin:Indxmax)=depthI;
 Depth_influx_line(1:Indxmin-1)=depthI(1);
 Depth_influx_line(Indxmax+1:end)=depthI(end);
 Yline=zeros(dom.Nx,1);   
 Yline(Indxmin:Indxmax)=XYinfl(:,2);
 Yline(1:Indxmin-1)=XYinfl(1,2);
 Yline(Indxmax+1:end)=XYinfl(end,2);
 Yline=repmat(Yline',Nomg,[]);
 
 [X,WW]=meshgrid(dom.X,halfomsig);
 KK=funOprt_invOmExact(dom.kx,halfomsig,Depth_influx_line,g);
 Rtheta=repmat(Rtheta',1,dom.Nx);
 RP=repmat(RP,1,dom.Nx);
 a_bar=repmat(a_bar',1,dom.Nx);
 
  for ii=1:length(timesig)          
      eta_x_kjj_thetakk=a_bar.*cos(KK.*sin(Rtheta+mDir_rad).*Yline...
                                +KK.*cos(Rtheta+mDir_rad).*X...
                                -WW.*timesig(ii)+RP);
      INsig(ii,:)=sum(eta_x_kjj_thetakk);
  end
end

