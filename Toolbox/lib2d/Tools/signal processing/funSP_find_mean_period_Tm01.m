function Tm01=funSP_find_mean_period_Tm01(T,Eta)
[oms,Sww]=funSP_variance_density_spectrum1d(T,Eta);
Tm01=2*pi*trapz(oms,Sww)/trapz(oms,oms.*Sww);