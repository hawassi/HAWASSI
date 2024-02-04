function tauG=fun_tangential_flow_global(k,F0,chiS,eta,Nonlin)
if Nonlin<=2
Nonlin=2;    
end
tauG=Ifft(funOprt_Finv(k,F0,eta,fft(chiS),Nonlin));