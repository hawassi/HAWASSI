function KK=funOprt_invOmBBM(k,ww,Depth,g)

 k_half        =k(1:floor(end/2));
 KK=zeros(length(ww),length(Depth));
 for ii=1:length(Depth) 
 Omega          = funOprt_OmBBM(k_half,Depth(ii),g,[]);
 KK(:,ii)       = interp1(Omega,k_half,ww,'spline');
 end 