function JS_NOR=funG_JONSWAP_1D_spectrum(g,gam,omg,Tp,Hs)


%%%% JONSWAP formula based on : Waves in Oceanic and coastal waters by Leo H. Holthuijsen
% This is non-normalized spectra !!
omgp                = 2*pi./Tp;
sigma               = zeros(size(omg));
sigma(omg<=omgp)    = 0.07;
sigma(omg> omgp)    = 0.09;
alpha               = 0.0081;

f           = omg./(2*pi);
f_p         = omgp./(2*pi);
gampangkat  = exp( -0.5.*( (f./f_p - 1)./sigma).^2 );
%JPM         = alpha.*g.^2.*((2*pi).^(-4)).*(f.^(-5)).*exp( -(5/4).*(f./f_p).^(-4)  ) ; %Pierson - Moskowitz
JS          =  alpha.*g.^2.*((2*pi).^(-4)).*(f.^(-5)).*exp( -(5/4).*(f./f_p).^(-4)  ) .* gam.^(gampangkat); %JONSWAP

% calculating Hs from non-normalized JONSWAP spectrum :
m0      = trapz(omg,JS);  % mo = int(omg,JS(omg));
%Hm0     = 4*sqrt(m0);             

JS_NOR  = (Hs^2/16).*JS./m0;
