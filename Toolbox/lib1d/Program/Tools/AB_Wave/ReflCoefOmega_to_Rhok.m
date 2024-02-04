function [Rho_k,Rho_Om,R_Om,R_k]=ReflCoefOmega_to_Rhok(k,D,Om,omAdd,ww,par)

%rho_Om=tanh(om);
%rho_Om=min(om.^2,1);
fun=['@(f)(',par.wall.file,')'];
R_Om_fun=str2func(fun);
R_Om=R_Om_fun(ww);


Rho_Om=(1-R_Om)./(1+R_Om);
k_Om=interp1(Om(k,D,omAdd),k,ww,'spline');
Rho_k=interp1(k_Om,Rho_Om,k,'spline');
R_k=interp1(k_Om,R_Om,k,'spline');

% figure
% plot(ww(1:floor(end/2)),refl_coef(1:floor(end/2)));
% xlabel('\omega(k)')
% ylabel('R')
% ylim([0 1])