function [Dmin,Dmax,gammin,gammax]= DispIP(k,Up,nupeak,bathy)
% 2-parameter approx of phase velocity by interpolation
% Capp(k,D) =gammin(D)*Up(k,Dmin)+gammax(D)*Up(k,Dmax)
% parpool(2);
% Up = str2func('UpExact');
        Bathy   = -bathy;
        Dmin    = min(Bathy);
        Dmax    = max(Bathy);
        N       = length(Bathy);
        Cmin_0  = sqrt(g*Dmin);
        Cmax_0  = sqrt(g*Dmax);
        kapnu   = zeros(N,1);
parfor j= 1:N
    kapnu(j)    = interp1(Om(k,Bathy(j)),k,nupeak); %kapnu(D) = k at peak frequency nu
end
        Cminnu  = Up(kapnu,Dmin)';
        Cmaxnu  = Up(kapnu,Dmax)';
        Cnuv    = zeros(N,1)';
parfor j= 1:N
    Cnuv(j)    = UpExact(kapnu(j),Bathy(j));%feval(Cnu,j);%interp1(Om(k,Bathy(j)),k,nupeak); %kapnu(D) = k at peak frequency nu
end
        det     = (Cmin_0.*Cmaxnu) -(Cmax_0.* Cminnu);
        gammin  = ( (Cmaxnu.* sqrt(g*Bathy))-(Cmax_0.*Cnuv) )./det;
        gammax = ( (Cmin_0.*Cnuv)-(Cminnu.*sqrt(g*Bathy)) )./det;      
figure(100)
plot(Bathy,gammin,'b', Bathy, gammax,'r')        
end
% err   =zeros(length(Bathy),length(k));
% errknu=zeros(length(Bathy),1);
% knu   =kapnu;%zeros(length(Bathy),1);
%  parfor j=1:length(Bathy)
%      for m=1:N
%          kkk=k(m);
%      Cp_exact=(Up(kkk,Bathy(j)));
%      Cp_approx=gammin(j)*(Up(kkk,Dmin))+...
%                 gammax(j)*(Up(kkk,Dmax)); 
%      err(j,:)=abs(Cp_exact-Cp_approx);
% %      knu(j)= interp1(Om(kkk,Bathy(j)),kkk,nupeak);
%      indknu  =closest(kkk,knu(j));
%      errknu(j)=abs(Cp_exact(indknu)-Cp_approx(indknu));
%      end;
%  end
%  
% figure
% mesh(Bathy,k,err');
% view(0,90);
% xlabel('depth (m) ');
% ylabel('wave number');
% axis([0 max(Bathy) 0 max(k)]);
% hold on;
% h=colorbar;
% %caxis([0 1]);
% title('Error of Interpolation of Cp')
% ylabel(h, 'Cp');
% lineobj = findobj('type', 'line');
% set(lineobj, 'linewidth', 1.3);
% textobj = findobj('type', 'text');
% set(textobj, 'fontunits', 'points');
% set(textobj, 'fontsize', 10);
% set(textobj, 'fontweight', 'bold');
% saveas(gcf,['ErrIn3p'],'fig');
% 
%  %Error at wavenumber in the peak frequency
% figure(10)
% plot(Bathy,errknu,'r',Bathy,Cnuv,'b');
% title('Error (red) and Cp at \kappa(x)')
% xlabel('depth (m)');
% ylabel('Cp');
% grid on;    
% lineobj = findobj('type', 'line');
% set(lineobj, 'linewidth', 1.3);
% textobj = findobj('type', 'text');
% set(textobj, 'fontunits', 'points');
% set(textobj, 'fontsize', 10);
% set(textobj, 'fontweight', 'bold');
% saveas(gcf,['errknu'],'fig');          
% 
% figure
% plot(Bathy,knu);
% title('wave number at peak frequency');
% xlabel('depth');
% ylabel('\kappa');
