function [R_Om,R_k]=funW_ReflCoefOmega_to_RhoK(g,dom,Om,omAdd,ww)

%rho_Om=tanh(om);
%rho_Om=min(om.^2,1);
for II=1:dom.wall.N
    if dom.wall.ReflCoef.FreqDep(II).flag==1
        fun=['@(w)(',cell2mat(dom.wall.param(II,12)),')'];
        R_Om_fun=str2func(fun);
        R_Om(II).coef=R_Om_fun(ww);
        
        [kkx,kky]=meshgrid(dom.kx,dom.ky);
        KK=sqrt(kkx.^2+kky.^2);
        K0=linspace(min(min(KK)),max(max(KK)),1000)';
        Depth=mean(mean(-dom.bathy.profile(dom.wall.ReflCoef.FreqDep(II).char==1)));
        %Rho_Om=(1-R_Om)./(1+R_Om);

        K_Om=interp1(Om(K0,Depth,g,omAdd),K0,ww,'spline');
        
        %Rho_k=interp1(k_Om,Rho_Om,k,'spline');
        R_k(II).coef=interp1(K_Om,R_Om(II).coef,KK,'spline');
    else
        R_Om(II).coef=0;  R_k(II).coef=0;
    end
end

% figure
% plot(ww(1:floor(end/2)),R_Om(1).coef(1:floor(end/2)));
% xlabel('\omega(k)')
% ylabel('R')
% ylim([0 1])