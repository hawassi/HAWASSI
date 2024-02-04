VarDens=load('target_spectrum_T01_2.txt');
ff=VarDens(2:end,1);
VarD=VarDens(2:end,2);
Hs=3;
M0=trapz(ff,VarD)*2*pi;
VarD=(Hs^2/16)*VarD./M0;
h0=4*sqrt(trapz(ff,VarD)*2*pi);

VarDensSpect=[2*pi*ff VarD];
save VarDensSpect VarDensSpect
figure(5)
plot(2*pi*ff,VarD)
xlabel('\omega [rad/s]');
ylabel('variance density [m^2 s/rad]');
