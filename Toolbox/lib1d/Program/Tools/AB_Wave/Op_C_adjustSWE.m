function Cvar_hat=Op_C_adjustSWE(g,H,H_min,Cvar_hat,var)
Cvar=Ifft(Cvar_hat);
Cvar(H<H_min)=sqrt(g*H(H<H_min)).*var(H<H_min);
Cvar_hat=fft(Cvar);