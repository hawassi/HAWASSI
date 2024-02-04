function Kx=localWaveNumber(eta,x)
Hx=imag(hilbert(eta)); dx=x(2)-x(1);
Kx=(eta.*gradient(Hx,dx)-Hx.*gradient(eta,dx))./(eta.^2+Hx.^2);