function [gam_min,gam_plus]= sinhcosh_2IP(D,D_min,D_plus,k,Om,nu)
% 2-parameter approx of phase velocity by interpolation
% (from Ruddy Beji-Battjes): cond: Om exact and C^2 exact at kapnu for disp of  C^2
% operator
g=9.81;
kappanu_D   = zeros(length(D),1);
cosh_D      = zeros(length(D),1);
cosh_min    = zeros(length(D),1);
cosh_plus   = zeros(length(D),1);
sinh_D      = zeros(length(D),1);
sinh_min    = zeros(length(D),1);
sinh_plus   = zeros(length(D),1);

for j= 1:length(D)
    if D(j)==0
    kappanu_D(j)      = k(2);    
    else
    kappanu_D(j)      = funOprt_invOmExactVal(nu,D(j),g);%interp1(Om(k,D(j)),k,nu);%interpolation to find k at peak frequency nu
    end
    cosh_min(j)       = cosh(kappanu_D(j)*D_min);
    cosh_plus(j)      = cosh(kappanu_D(j)*D_plus);
    cosh_D(j)         = cosh(kappanu_D(j)*D(j));
    sinh_min(j)       = sinh(kappanu_D(j)*D_min);
    sinh_plus(j)      = sinh(kappanu_D(j)*D_plus);
    sinh_D(j)         = sinh(kappanu_D(j)*D(j));
end

    det         = cosh_min.*sinh_plus -cosh_plus.*sinh_min;
    gam_min     = ( sinh_plus.*cosh_D-cosh_plus.*sinh_D )./det;
    gam_plus    = ( -sinh_min.*cosh_D+cosh_min.*sinh_D)./det;
    
%  assignin('base','D',D) ;   
% hf1=figure;
% plot(D,gam_min,'b', D, gam_plus,'r')  
% set(hf1,'units','normalized','Position',[0.1 0.1 0.5 0.5])
% pause;

