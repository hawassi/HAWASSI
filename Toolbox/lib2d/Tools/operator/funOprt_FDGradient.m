function [Gradx,Grady]=funOprt_FDGradient(M,dx,dy,Id)
if Id==1 % Forward difference
 Mxplus1=circshift(M,-1,2);
 Myplus1=circshift(M,-1,1);
 Gradx=(Mxplus1-M)./dx;
 Grady=(Myplus1-M)./dy;
 Gradx(:,1)=0;Gradx(:,end)=0;
 Grady(1,:)=0;Grady(end,:)=0;
elseif Id==-1 % Backward difference
 Mxmin1=circshift(M,1,2);
 Mymin1=circshift(M,1,1);
 Gradx=(M-Mxmin1)./dx;
 Grady=(M-Mymin1)./dy;   
 Gradx(:,1)=0;Gradx(:,end)=0;
 Grady(1,:)=0;Grady(end,:)=0;
else % central difference
 [Gradx,Grady]=gradient(M,dx,dy);   
end

