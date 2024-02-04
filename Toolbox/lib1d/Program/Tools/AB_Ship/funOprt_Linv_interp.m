function LinvVar_hat=funOprt_Linv_interp(g,kk,Oprt,var,bathtype)

Integvar_hat=fft(var)./(1i.*kk);
Integvar_hat(isnan(Integvar_hat))=0;
Integvar_hat(1)=0;
Integvar=Ifft(Integvar_hat);

LinvVar_hat=funOprt_Finv_Interp(g,kk,Oprt,-Integvar,bathtype);



% if OrdeOpL~=1
% L0invphi_hat=Linvphi_hat;   
% L0L0invphi_hat=L0.*L0invphi_hat;
% L0L0invphi_hat(isnan(L0L0invphi_hat))=0;
% dxL0invphi_hat=1i.*kk.*L0invphi_hat;
% dxL0invphi_hat(isnan(dxL0invphi_hat))=0;
% M1L0invphi_hat=-L0.*fft(zeta.*Ifft(L0L0invphi_hat))-1i.*kk.*fft(zeta.*Ifft(dxL0invphi_hat));
% term2order=L0inv.*M1L0invphi_hat;
% Linvphi_hat=Linvphi_hat+term2order;
% if OrdeOpL~=2
% term3order=-0.5.*(L0inv.*(1i.*kk).^2.*fft(zeta.^2.*var)+fft(zeta.^2.*Ifft((1i.*kk).^2.*L0inv.*varhat)))...
%     +fft(zeta.*Ifft(1i.*kk.*fft(zeta.*Ifft(1i.*kk.*L0inv.*varhat))))...
%     +L0inv.*1i.*kk.*fft(zeta.*Ifft(1i.*kk.*fft(zeta.*var)))+...
%     +L0inv.*1i.*kk.*fft(zeta.*Ifft(1i.*kk.*L0inv.*1i.*kk.*fft(zeta.*Ifft(1i.*kk.*L0inv.*varhat))));
% Linvphi_hat=Linvphi_hat+term3order;
% end
% end
% Linvphi_hat(1)=0;
% end