function dynmodel = DynModelName(model,bath)
evolution     =model.evol;
parnonlinear  =model.nonlinear;
breaking      =model.breaking.check;
bathtype      =bath.type;
wallpresence  =bath.wall.check;

if    strcmp(breaking,'No')     
    if strcmp(wallpresence,'Yes')
        dynmodel   = [evolution,num2str(parnonlinear),bathtype,'Wall'];
    else
        dynmodel   = [evolution,num2str(parnonlinear),bathtype]; 
    end;
else %    strcmp(breaking,'Yes')     
    if strcmp(wallpresence,'Yes')
        dynmodel   = [evolution,num2str(parnonlinear),'br',bathtype,'Wall'];
    else
        dynmodel   = [evolution,num2str(parnonlinear),'br',bathtype]; 
    end
end
