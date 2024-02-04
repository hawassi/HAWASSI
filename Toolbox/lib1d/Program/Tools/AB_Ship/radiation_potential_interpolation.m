function radpar=radiation_potential_interpolation(shippar,par,model,influx,Proj)

chiship=shippar.form.chiXcZ0(:,end);
zeta=shippar.form.shapeXcZ0(end,:)';

Nship=shippar.Nship;
nu=0;
for ii=1:Nship
    if strcmp(shippar.data(ii,2),'Heave')
    nu=nu+shippar.form.nuXcZ0.z(ii,:)';
    elseif strcmp(shippar.data(ii,2),'Surge')
    nu=nu+shippar.form.nuXcZ0.x(ii,:)';  
    elseif strcmp(shippar.data(ii,2),'Fixed')
    chiship=chiship.*(1-shippar.form.chi(:,ii));
    end
end


 
  radpar=[];

end