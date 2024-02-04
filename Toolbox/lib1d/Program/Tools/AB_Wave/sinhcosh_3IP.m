function [gam_min,gam_mid,gam_plus]= sinhcosh_3IP(Htot,H_min,H_plus,H_mid,k,Om,nu)
Ninterp=length(Htot);
Cosh_H = zeros(Ninterp,1);
Cosh_min = zeros(Ninterp,1);
Cosh_plus = zeros(Ninterp,1);
Cosh_mid = zeros(Ninterp,1);
Sinh_H = zeros(Ninterp,1);
Sinh_min = zeros(Ninterp,1);
Sinh_plus = zeros(Ninterp,1);
Sinh_mid = zeros(Ninterp,1);
Tanh_H = zeros(Ninterp,1);
Tanh_min = zeros(Ninterp,1);
Tanh_plus = zeros(Ninterp,1);
Tanh_mid = zeros(Ninterp,1);


for j= 1:Ninterp
kappanu_H   = funOprt_invOmExactVal(nu,Htot(j),g);%interp1(Om(k,Htot(j)),k,nu);%interpolation to find k at peak frequency nu
Cosh_min(j) = cosh(kappanu_H.*H_min);
Cosh_plus(j)= cosh(kappanu_H.*H_plus);
Cosh_mid(j) = cosh(kappanu_H.*H_mid);
Cosh_H(j)   = cosh(kappanu_H.*Htot(j));

Sinh_min(j) = sinh(kappanu_H.*H_min);
Sinh_plus(j)= sinh(kappanu_H.*H_plus);
Sinh_mid(j) = sinh(kappanu_H.*H_mid);
Sinh_H(j)   = sinh(kappanu_H.*Htot(j));

Tanh_min(j) = tanh(kappanu_H.*H_min);
Tanh_plus(j)= tanh(kappanu_H.*H_plus);
Tanh_mid(j) = tanh(kappanu_H.*H_mid);
Tanh_H(j)   = tanh(kappanu_H.*Htot(j));
end


A= (Sinh_mid.*Tanh_plus)-(Sinh_plus.*Tanh_mid);
B=-(Sinh_min.*Tanh_plus)+(Sinh_plus.*Tanh_min);
C= (Sinh_min.*Tanh_mid)-(Sinh_mid.*Tanh_min);
D=-(Cosh_mid.*Tanh_plus)+(Cosh_plus.*Tanh_mid);
E= (Cosh_min.*Tanh_plus)-(Cosh_plus.*Tanh_min);

F=-(Cosh_min.*Tanh_mid)+(Cosh_mid.*Tanh_min);
G= (Cosh_mid.*Sinh_plus)-(Cosh_plus.*Sinh_mid);
H=-(Cosh_min.*Sinh_plus)+(Cosh_plus.*Sinh_min);
I= (Cosh_min.*Sinh_mid)-(Cosh_mid.*Sinh_min);


det=Cosh_min.*A+Cosh_mid.*B+Cosh_plus.*C;

gam_min  = ( A.*Cosh_H+D.*Sinh_H+G.*Tanh_H )./det;
gam_mid  = ( B.*Cosh_H+E.*Sinh_H+H.*Tanh_H )./det;
gam_plus = ( C.*Cosh_H+F.*Sinh_H+I.*Tanh_H )./det;
    
    
    
    
    
   
% hf1=figure;
% plot(Htot,gam_min,'b',Htot, gam_plus,'r',Htot, gam_mid,'g')  
% set(hf1,'units','normalized','Position',[0.1 0.1 0.5 0.5])
% pause;

