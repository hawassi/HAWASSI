function plot_PP1_V(X,Z,V,Vx,Vz,string,labels,CheckColorbar,Clim)
%surfl(X,Z,Var);colormap('winter'); shading interp
mesh(X,Z,V);hold on;
quiver(X,Z,Vx,Vz,'color',[0 0 1]);
xlim([min(X(1),X(end)) max(X(1),X(end))]);
ylim([min(Z(1),Z(end)) max(Z(1),Z(end))]);
xlabel(labels{1});
ylabel(labels{2});
if CheckColorbar==1
caxis([Clim(1) Clim(2)])
end
view(2);
h=colorbar;
ylabel(h,string)