function tauG=fun_tangential_flow_global_Interp(g,k,Oprt,varS,bathtype)

tauG=Ifft(funOprt_Finv_Interp(g,k,Oprt,varS,bathtype));