function taunorm=normalized_tangential_potential_interpolation(shippar,par,model)

chiship=shippar.form.chiXcZ0(:,end);
chiship_hat=fft(chiship);
Cp2inv_minChiship_hat=shippar.Oprt.Cp2invmin.*chiship_hat;
Cp2inv_minChiship=Ifft(Cp2inv_minChiship_hat);
Cp2inv_plusChiship_hat=shippar.Oprt.Cp2invplus.*chiship_hat;
Cp2inv_plusChiship=Ifft(Cp2inv_plusChiship_hat);

taunorm=HSS(A,a,chiship_hat,chiship)+HSS();

end