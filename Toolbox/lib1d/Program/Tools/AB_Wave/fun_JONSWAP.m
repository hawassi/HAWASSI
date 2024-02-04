function [JS_out,omg,omgp,domg] = fun_JONSWAP(Nomg,lT,Tp,Hs,gam,g,plotting)
% ========================================================================%
% function to create 1D JONSWAP spectra
% The calculation is based on "Waves in Oceanic and coastal waters" by 
% Leo H. Holthuijsen
% The spectra is normalized by using the definition of zeroth order moment
% m0 = int(omg,JONSWAP) and its relation with significant wave height Hs :
% m0 = Hs^2/16, or Hs = 4*sqrt(m0);
% 
% Copyright (c) LabMath-Indonesia & University of Twente 
%          didit@labmath-indonesia.org 
%          didit.adytia@gmail.com
%          v.02. 23 May 2013
% Remark :
% - The normalization factor for Hs is NOW correct yet for 2D-spectra, use
% the definition of zeroth order moment : mo = Hs^2/16.
% ========================================================================%

domg   = 2*pi/lT;                     %delta omega
omg    = [domg : domg : (Nomg)*domg ];   %omega vector with Nomg # of discretization in omega 
omgp   = 2*pi/Tp;                        %omega peak

%% JONSWAP formula based on : Waves in Oceanic and coastal waters by Leo H. Holthuijsen
% notice that this is non-normalized spectra !!
sigma               = zeros(size(omg));
sigma(omg<=omgp)    = 0.07;
sigma(omg> omgp)    = 0.09;
a                   = exp( - (omg - omgp).^2./(2.*omgp.^2*sigma.^2 )   );
beta                = -5/4;
alpha               = 0.0081;

JS          = zeros(size(omg));
f           = omg./(2*pi);
f_p         = omgp./(2*pi);
gampangkat  = exp( -0.5.*( (f./f_p - 1)./sigma).^2 );
JPM         = alpha.*g.^2.*((2*pi).^(-4)).*(f.^(-5)).*exp( -(5/4).*(f./f_p).^(-4)  ) ; %Pierson - Moskowitz
JS          =  alpha.*g.^2.*((2*pi).^(-4)).*(f.^(-5)).*exp( -(5/4).*(f./f_p).^(-4)  ) .* gam.^(gampangkat); %JONSWAP

% calculating Hs from non-normalized JONSWAP spectrum :
m0      = trapz(omg,JS);  % mo = int(omg,JS(omg));
Hm0     = 4*sqrt(m0);             
Hs_v1   = 1*Hm0;

% Now normalize the spectrum based on definition of zeroth-order moment 
% and its relation with Hs : mo = int(omg,JS(omega)) and Hs = 4*sqrt(mo); so:

JS_NOR  = (Hs^2/16).*JS./(trapz(omg,JS));

if strcmp(plotting,'Yes')
%     figure;
%     subplot(2,1,1);
%     plot(omg,JS); xlabel('\omega [Rad/s]'); ylabel('JONSWAP');
%     axis([0 , 8*omgp, 0 , 1.1*max(max(JS))]);
%     title('JONSWAP before normalization with Hs');
%     subplot(2,1,2);
%     plot(f,JS); xlabel('f [Hz]'); ylabel('JONSWAP')
%     axis([0 , 8*f_p, 0 , 1.1*max(max(JS))])

    figure
        subplot(2,1,1);
        plot(omg,JS_NOR); xlabel('\omega [Rad/s]'); ylabel('JONSWAP');
        title(['JONSWAP spectrum, Hs= ',num2str(Hs),', Tp = ',num2str(Tp)]);
        subplot(2,1,2);
        plot(f,JS_NOR); xlabel('f [Hz]'); ylabel('JONSWAP')
end
% calculate Hs :
m0      = trapz(omg,JS_NOR);
Hm0     = 4*sqrt(m0);
Hs_v2   = 1*Hm0;
% saving output :
JS_out = JS_NOR;


end