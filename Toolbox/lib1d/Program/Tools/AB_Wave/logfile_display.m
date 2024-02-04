fc=fix(clock);
par.fc=fc;

% prev_log=false;
% if exist([Proj.Dir,Proj.savename,'.log'],'file')
% fid = fopen([Proj.Dir,Proj.savename,'.log'],'r');
% textLogPrev = textscan(fid,'%s','delimiter','\n','CollectOutput', true);
% assignin('base','text',textLogPrev)
% fclose(fid);
% delete([Proj.Dir,Proj.savename,'.log']);
% prev_log=true;
% end


diary([Proj.Dir,Proj.savename,'.log']);
diary on;
disp('*********************************************************************')
disp(['HAWASSI-AB simulation diary, date: ', num2str(par.fc(3)),'.' num2str(par.fc(2)),'.' num2str(par.fc(1)),...
    ' time:  ' num2str(par.fc(4)),'h.' num2str(par.fc(5)),'min.'])
disp('*********************************************************************')
disp('%');
disp(['Project name: ',Proj.savename]);
disp('%');
%disp('-------------------------- UserNote reads -------------------------: ');
disp(['User note   :',Proj.UseNote]);
%disp('-------------------------- End UserNote ---------------------------: ');
disp('%');
disp('MODEL DESCRIPTION :');
disp(['Dynamic Model        : ', model.dyn(1:3)]);
disp(['Dispersion Model     : ', model.dispersion(3:end)]);
if strcmp(model.breaking.check,'Yes')
    disp(['Breaking             : ']);
    disp(['   Initiation        : U/C= ', num2str(model.breaking.KBC)]);
    disp(['   Termination       : uF/uI= ', num2str(model.breaking.TC),', T*/Tp= ', num2str(model.breaking.Tchar)]);
else
    disp(['Breaking             : No']);
end
disp('%');

disp('INFLUX DESCRIPTION : ');
if input.type==1
    disp( ['Signal type         : ',influx.name]);
end
if input.type~=1
    disp(['Depth (h)                     : ', num2str(influx.depth),'[m]']);
    if strcmp(influx.name,'Harmonic')||strcmp(influx.name,'User-defined')
        disp(['Signal type                   : ',influx.name]);
    elseif strcmp(influx.name,'Jonswap')
        disp(['Signal type                   : ',influx.name,' with gamma= ',num2str(influx.JS_gamma)]);
    end
    if strcmp(influx.name,'Harmonic')
        disp(['Amplitude (A)                 : ',num2str(roundn(influx.Hs/2,-3)),'[m]']);
    else
        disp(['Significant wave Height (Hs)  : ',num2str(roundn(influx.Hs,-3)),'[m]']);
    end
    disp(['Peak period (Tp)              : ',num2str(roundn(influx.Tp,-3)),'[s]']);
    disp(['      Derived info:']);
    disp(['      Peak frequency (nu)           : ',num2str(roundn(influx.nu_p,-3)),'[rad/s]']);
    disp(['      Peak wave-number (kp)         : ',num2str(roundn(influx.k_p,-3))]);
    disp(['      Peak wave-length              : ',num2str(roundn(influx.lambda_p,-3)),'[m]']);
    disp(['      Peak phase speed              : ',num2str(roundn(influx.Cp_p,-3)),'[m/s]']);
    disp(['      Peak group speed              : ',num2str(roundn(influx.Vg_p,-3)),'[m/s]']);
    disp(['      Steepness (kp*(Hs./2))        : ',num2str(roundn(influx.ka,-3))]);
    disp(['      Relative wave-length(lambda/h): ', num2str(influx.lambda_per_H),'(',influx.category,')']);
    disp(['                          (kp*h)    : ',num2str(roundn(influx.kh,-3))]);
end
disp('%');
disp('INITIAL WAVE CONDITIONS : ');
disp( ['Initial condition     : ',IVP.typename]);
if  IVP.type==2||IVP.type==3
    disp( ['Amplitude             : ',num2str(IVP.A)]);
    disp( ['Standard deviation    : ',num2str(IVP.lambda)]);
    disp( ['Center position       : ',num2str(IVP.x0)]);
end
if IVP.type==4
    disp( ['File name             : ',IVP.filename]);
end
disp('%')
if influx.wind.check==1
    disp('WIND INPUT : YES');
    disp(['     Coef  : ', num2str(influx.wind.coef)]);
else
    disp('WIND INPUT : NO');
end
disp('%');
disp('NUMERICAL SETTINGS : ');
disp(['Spatial interval      : ','(', num2str(par.Xleft),';',num2str(par.Xright),') [m]']);
disp(['Length of Fourier bdy [L;R] : [',num2str(par.bf0(1)),';',num2str(par.bf0(2)),'] [m]']);
disp(['Number of Nodes       : ',num2str(par.Nx)]);
disp(['        Grid size (dx): ',num2str(roundn(par.dx,-3)),'[m]']);
disp(['        Cutfrac k     : ',num2str(par.cutfrac)]);
disp(['Time interval         : ','(', num2str(par.t_init) ,';',num2str(par.t_end),') [s]']);
disp(['        Time step (dt): ', num2str(roundn(par.dt,-3)),' [s]']);
disp('%');
if strcmp(bath.type,'F')
    disp(['Bathymetry            : ', bath.name]);
    disp(['     Depth            : ', num2str(bath.depth),'[m]']);
elseif strcmp(bath.type,'B')
    disp(['Bathymetry            : ', bath.name])
    disp(['     Depth            : max= ', num2str(bath.par(1)),'[m],',' min= ', num2str(bath.par(2)),'[m]']);
    disp(['     Slope            : gradient= ', num2str(bath.par(3)), ', start position= ', num2str(bath.par(4)),'[m]']);
elseif strcmp(bath.type,'BR')
    disp(['Bathymetry            : ', bath.name])
    disp(['     Depth            : max= ', num2str(bath.par(1)),'[m]']);
    disp(['     Slope            : gradient=', num2str(bath.par(2)), ', shore position : ', num2str(bath.par(3)),'[m]']);
else
    disp(['Bathymetry            : ', bath.name]);
    disp(['     File-name        : ', bath.filename]);
end
disp('%');

if bath.friction.check==1
    strF=bath.friction.data;
    NsF=length(strF.cf);
    for i=1:NsF-1
        strCf{i}=[strF.cf{i},' , ',];
        strInterv{i}=['(',strF.interval{i},'),'];
    end
    strCf{NsF}=[strF.cf{NsF}];
    strInterv{NsF}=['(',strF.interval{NsF},')'];
    disp(['Bottom friction: ']);
    disp(['      coef         : ', strCf{1:NsF}]);
    disp(['      interval     : ',strInterv{1:NsF},' [m]']);
else
    disp(['Bottom friction: No']);
end

disp('%');
if strcmp(bath.wall.check,'Yes')
    if bath.wall.type==1
        disp(['Wall: ']);
        disp(['    Type       : Uniform']);
        disp(['    Position   : ',num2str(bath.wall.position), '[m]']);
        disp(['    Coef       : ',num2str(bath.wall.Coef)]);
    else
        disp(['Wall:']);
        disp(['    Type       : Frequency dependent']);
        disp(['    Position   : ',num2str(bath.wall.position), '[m]']);
        disp(['    Function  : R(f)=', bath.wall.file_R]);
    end
else
    disp(['Wall: No']);
end
disp('%');

if ~strcmp(model.influx.direction,'None')
    if strcmp(model.influx.type,'Point')
        disp(['Influx generation method : ', model.influx.direction,'; ', model.influx.type]);
    else
        disp(['Influx generation method : ', model.influx.direction,'; ', model.influx.type(3:end)]);
    end
    
    if ~isempty(bath.influx_AdjZone)
        disp(['    Nonlinear adjustment : ', num2str(bath.influx_AdjZone), '*lambda_peak']);
    else
        disp(['    Nonlinear adjustment : 0*lambda_peak']);
    end
    if input.ramp.check==1
        disp(['    Ramp                 : ', num2str(input.ramp.length), '*Tp']);
    else
        disp(['    Ramp                 : Off']);
    end
    if input.filter.check==1
        disp(['    Filter               : On; Freq. interval: (',num2str(input.filter.LFreq),';',num2str(input.filter.HFreq),') rad/s']);
    else
        disp(['    Filter               : Off']);
    end
    
    disp(['    Influx position      : ', num2str(influx.position),'[m]']);
else
    disp(['Influx generation method : ', model.influx.direction]);
end
disp('%');

if shippar.check==1
    
    disp('SHIP SETUP ');
    disp(['Nship: ',num2str(shippar.Nship)]);
    for ii=1:shippar.Nship
        if strcmpi(shippar.data(ii,1),'Barge')
            shipshape='Barge';
        elseif strcmpi(shippar.data(ii,1),'Half-circle')
            shipshape='Half-circle';
        else
            shipshape='User-defined';
        end
        
        if strcmp(shippar.data(ii,2),'Fixed')
            shipMot='Fixed';
            AdMass=0;
        elseif strcmp(shippar.data(ii,2),'Heave')
            shipMot='Heave';
            AdMass=shippar.rad.Ma.z(ii);
        elseif strcmp(shippar.data(ii,2),'Surge')
            shipMot='Surge';
            AdMass=shippar.rad.Ma.x(ii);
        elseif strcmp(shippar.data(ii,2),'Pitch')
            shipMot='Pitch';
            AdMass=shippar.rad.Ma.theta(ii);
        else
            shipMot='Free';
            AdMassx=shippar.rad.Ma.x(ii);
            AdMassz=shippar.rad.Ma.z(ii);
            AdMasstheta=shippar.rad.Ma.theta(ii);
        end
        
        disp(['Ship #',num2str(ii),': Shape=',shipshape, ', Motion=',shipMot,...
            ', Length=',num2str(cell2mat(shippar.data(ii,3))),'[m] , Draft=',num2str(cell2mat(shippar.data(ii,4))),'[m]']);
        disp(['      Initial center ship position(Xc,Zc)= (',num2str(cell2mat(shippar.data(ii,5))),',',num2str(cell2mat(shippar.data(ii,6))),')[m]']);
        if strcmpi(shipMot,'Free')
            disp(['      Mass=', num2str(shippar.form.Mass(ii)), 'Moment Inertia=',num2str(shippar.form.MIner(ii)),...
                ', Added Mass (heave)= ', num2str(AdMassz),'(surge)= ', num2str(AdMassx),'(pitch)= ', num2str(AdMasstheta)]);
        elseif strcmpi(shipMot,'Pitch') 
            disp(['      Moment Inertia=', num2str(shippar.form.MIner(ii)),...
                ', Added Moment Inertia = ', num2str(AdMass)]);
        else
            disp(['      Mass=', num2str(shippar.form.Mass(ii)),...
                ', Added Mass = ', num2str(AdMass)]);
        end
    end
else
    disp('SHIP:  No');
end
disp('%');
if options.interior.check==1
    disp(['INTERIOR CALCULATION PREPARATION: ']);
    disp(['Time interval        : (',num2str(options.interior.time(1)),';',num2str(options.interior.time(2)),') [s]']);
    disp(['Time step            : ',num2str(par.dt.*options.interior.time(3)),' [s]']);
    disp('%');
end

disp('ODE SETTING:');
disp(['ODE solver           : ', par.odesolver]);
disp(['Relative error tol.  : ', num2str(par.odetol)]);
disp(['No of Partition      : ', num2str(par.ode_Npartition)]);
disp('---------------------------------------------------------------------');



%%diary off; in HaWaSSI.m in function of numerical_setting_advises

%     if prev_log
%     disp(textLogPrev{1,1});
%     end
% if ~isdeployed
% clc;
% end
