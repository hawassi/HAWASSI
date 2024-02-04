function C2var_hat=Op_C2_adjustSWE(g,H,H_min,C2var_hat,var)
C2var=Ifft(C2var_hat);
C2var(H<H_min)=(g*H(H<H_min)).*var(H<H_min);
C2var_hat=fft(C2var);