function KK=funOprt_invOmExact(k,ww,Depth,g)
if abs(max(Depth)-min(Depth))<1e-3
   KK       =  funOprt_invOmExactVal(ww,max(Depth),g);
   KK       =  repmat(KK,1,length(Depth));
else
  %  KK      =  funOprt_invOmExactVal(ww,Depth,g);
    KK    =zeros(length(ww),length(Depth));
    ww_max=max(ww);
    for ii=1:length(Depth)
        kk_max=funOprt_invOmExactVal( ww_max,Depth(ii),g);
        kk=linspace(-1.1*kk_max,1.1*kk_max,10000);
        Omega          = funOprt_OmExact(kk,Depth(ii),g,[]);
        KK(:,ii)       = interp1(Omega,kk,ww,'spline');   %
    end
    
end

