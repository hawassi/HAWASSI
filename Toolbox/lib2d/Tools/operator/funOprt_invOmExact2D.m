function KK=funOprt_invOmExact2D(e,kin,ww,Depth,g) 

%%% still not finished
 %k_half        =k(1:floor(end/2));
 KK=zeros(length(ww),length(Depth));
 ww_max=max(ww);
 
 for ii=1:length(Depth) 
  kk_max=funOprt_invOmExactVal( ww_max,Depth(ii),g);
  kk=linspace(0,1.1*kk_max,200);
  
  Omega          = funOprt_Om2dExact(kk,kin,depth,g,e);%funOprt_OmExact(kk,Depth(ii),g,[]);
  
  KK(:,ii)       = interp1(Omega,kk,ww,'spline');   % 
 end
 