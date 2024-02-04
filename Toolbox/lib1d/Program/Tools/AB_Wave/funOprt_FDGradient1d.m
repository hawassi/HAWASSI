function [Gradx]=funOprt_FDGradient1d(M,dx,Id)
if Id==1 % Forward difference
 Mxplus1=circshift(M,-1);
 Gradx=(Mxplus1-M)./dx;
 Gradx(1)=0;Gradx(end)=0;
elseif Id==-1 % Backward difference
 Mxmin1=circshift(M,1);
 Gradx=(M-Mxmin1)./dx; 
 Gradx(1)=0;Gradx(end)=0;
else % central difference
 [Gradx]=gradient(M,dx);   
end

