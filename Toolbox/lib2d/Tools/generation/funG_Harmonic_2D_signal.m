function [INsig,SpectProp]=funG_Harmonic_2D_signal(g,ampl,Tp,dom,timesig,depthinfluxline,spatialinflux,theta0,I)

W0=2*pi./Tp;
%Depth=depthinfluxline;
XYinfl=spatialinflux.line(I).xy;
if strcmpi(spatialinflux.line(I).Orientation,'Horizontal')
    INsig=zeros(length(timesig),dom.Nx);
    
     infx1=funC_closest(dom.X,XYinfl(1,1));
     infxN=funC_closest(dom.X,XYinfl(end,1));
    
    for jj=1:dom.Nx%length(depthinfluxline)
        indxInfl=jj;
        if jj<infx1
        indyInfl= funC_closest(dom.Y,XYinfl(1,2));
        elseif jj>infxN
        indyInfl= funC_closest(dom.Y,XYinfl(end,2));
        else
         indyInfl=funC_closest(dom.Y,XYinfl(jj-infx1+1,2));
        end
%         Depthjj=-dom.bathy.profile(indyInfl,indxInfl);
%         K0=funOprt_invOmExact(dom.ky,W0,Depthjj,g);
%         
%         INsig(:,indxInfl)=ampl.*cos(K0.*cos(theta0)*dom.X(indxInfl)+K0.*sin(theta0)*dom.Y(indyInfl)-W0.*timesig);
  % nunu
        if jj == 1; yinit = dom.Y(indyInfl);  end	
        Depthjj = -dom.bathy.profile(indyInfl,indxInfl);
        K0      = funOprt_invOmExact(dom.ky,W0,Depthjj,g);
        INsig(:,indxInfl) = ampl.*cos(K0.*cos(theta0)*dom.X(indxInfl)+K0.*sin(theta0)*(dom.Y(indyInfl)-yinit)-W0.*timesig);
    end
else
    INsig=zeros(length(timesig),dom.Ny);
    
    infy1=funC_closest(dom.Y,XYinfl(1,2));
    infyN=funC_closest(dom.Y,XYinfl(end,2));
    
    for jj=1:dom.Ny %length(depthinfluxline)
        indyInfl=jj;
        if jj<infy1
        indxInfl= funC_closest(dom.X,XYinfl(1,1));
        elseif jj>infyN
        indxInfl= funC_closest(dom.X,XYinfl(end,1));
        else
        indxInfl=funC_closest(dom.X,XYinfl(jj-infy1+1,1));
        end
%         Depthjj=-dom.bathy.profile(indyInfl,indxInfl);
%         K0=funOprt_invOmExact(dom.kx,W0,Depthjj,g);
%         INsig(:,indyInfl)=ampl.*cos(K0.*cos(theta0)*dom.X(indxInfl)+K0.*sin(theta0)*dom.Y(indyInfl)-W0.*timesig);
 % nunu
		if jj == 1; xinit = dom.X(indxInfl); end
        Depthjj = -dom.bathy.profile(indyInfl,indxInfl);
        K0      = funOprt_invOmExact(dom.kx,W0,Depthjj,g);
        INsig(:,indyInfl)=ampl.*cos(K0.*cos(theta0)*(dom.X(indxInfl)-xinit)+K0.*sin(theta0)*dom.Y(indyInfl)-W0.*timesig);
    end
end
%  [YY,TT]=meshgrid(dom.Y,timesig);
%     figure;
%     surf(TT,YY,INsig,'edgecolor','none')
%     view(2)

%% for plotting 2d spectrum
omsig=funC_freqspace(timesig);
halfomsig=omsig(1:floor(end/2));
Nomg=length(halfomsig);                % uniform for long-crested waves
theta  = linspace(-pi/2,pi/2, Nomg );
StdTheta = 0;cos2pdf=zeros(size(halfomsig));
cos2pdf(funC_closest(theta,0))=1/(theta(2)-theta(1));

       if strcmp(spatialinflux.line(I).Orientation,'Vertical')
           Ymid=(max(XYinfl(:,2))+min(XYinfl(:,2)))/2;
           IndS=funC_closest(dom.Y,Ymid);
       else
           Xmid=(max(XYinfl(:,1))+min(XYinfl(:,1)))/2;
           IndS=funC_closest(dom.X,Xmid);
       end
       
[OmS,S_hat]=funSP_variance_density_spectrum1d(timesig,INsig(:,IndS));

SpectProp.varDensSpect=S_hat';
SpectProp.StdTheta=StdTheta;
SpectProp.cos2pdf=cos2pdf';
SpectProp.mDir=theta0;
SpectProp.halfomsig=OmS;