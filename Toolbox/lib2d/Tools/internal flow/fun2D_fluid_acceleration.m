function [dt_dxPhi,dt_dyPhi,dt_dzPhi]=fun2D_fluid_acceleration(dxPhi,dyPhi,dzPhi,T,Z,Y,X)
% dt_dxPhi=zeros(size(dxPhi));
% dt_dyPhi=zeros(size(dyPhi));
% dt_dzPhi=zeros(size(dzPhi));
% 
% for ii=1:length(Z)
%     for jj=1:length(Y)
%         for kk=1:length(X)
%             dxPhi_t=squeeze(squeeze(squeeze(dxPhi(:,ii,jj,kk))));
%             dyPhi_t=squeeze(squeeze(squeeze(dyPhi(:,ii,jj,kk))));
%             dzPhi_t=squeeze(squeeze(squeeze(dzPhi(:,ii,jj,kk))));
%             dt_dxPhi(:,ii,jj,kk)=gradient(dxPhi_t,T(2)-T(1));
%             dt_dyPhi(:,ii,jj,kk)=gradient(dyPhi_t,T(2)-T(1));
%             dt_dzPhi(:,ii,jj,kk)=gradient(dzPhi_t,T(2)-T(1));
%         end
%     end
%  end

[~,dt_dxPhi,~,~]=gradient(dxPhi,Z(2)-Z(1),T(2)-T(1),Y(2)-Y(1),X(2)-X(1));
[~,dt_dyPhi,~,~]=gradient(dyPhi,Z(2)-Z(1),T(2)-T(1),Y(2)-Y(1),X(2)-X(1));
[~,dt_dzPhi,~,~]=gradient(dzPhi,Z(2)-Z(1),T(2)-T(1),Y(2)-Y(1),X(2)-X(1));