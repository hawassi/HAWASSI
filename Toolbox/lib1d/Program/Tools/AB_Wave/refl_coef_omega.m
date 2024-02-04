function [ww,refl_coef,rho_Om]=refl_coef_omega(k,D,Om,par)

%rho_Om=tanh(om);
%rho_Om=min(om.^2,1);
ww=Om(k,D);
if par.wall.file_def==1
rho_Om=min(ww.^2,0.5);
else
rho_Om_fun=str2func(par.wall.file);
rho_Om=rho_Om_fun(ww);
end
Om_tilde=ww./rho_Om;
ktilde=invOm(Om_tilde,D,Om);
refl_coef=(ktilde-k)./(ktilde+k);
% figure
% plot(ww(1:floor(end/2)),refl_coef(1:floor(end/2)));
% xlabel('\omega(k)')
% ylabel('R')
% ylim([0 1])