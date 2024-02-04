function plot_PP1_Var(X,Z,Var,string,labels,SetAx,Id)
%surfl(X,Z,Var);colormap('winter'); shading interp


if SetAx.Xlim.Id==1
indx1=closest(X,SetAx.Xlim.val(1));
indx2=closest(X,SetAx.Xlim.val(2));
if Id==1 || Id==3
mesh(X(indx1:indx2),Z,Var(:,(indx1:indx2),:));
xlim([min(X(indx1),X(indx2)) max(X(indx1),X(indx2))]);
elseif Id==2
mesh(X(indx1:indx2),Z,Var((indx1:indx2),:,:));
xlim([min(X(indx1),X(indx2)) max(X(indx1),X(indx2))]);
end
else
mesh(X,Z,Var); 
xlim([min(X(1),X(end)) max(X(1),X(end))]);
end
if SetAx.Ylim.Id==1
indz1=closest(Z,SetAx.Ylim.val(1));
indz2=closest(Z,SetAx.Ylim.val(2));    
ylim([Z(indz1) Z(indz2)]);
else
ylim([min(Z(1),Z(end)) max(Z(1),Z(end))]);
end

xlabel(labels{1});
ylabel(labels{2});
if SetAx.Clim.Id==1
caxis([SetAx.Clim.val(1) SetAx.Clim.val(2)])
end
view(2);
h=colorbar;
ylabel(h,string)