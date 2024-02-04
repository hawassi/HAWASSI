% !--------------------------------------------------------------------------------------
% !
% !    Copyright (C) 2024 - LabMath-Indonesia
% !
% !    This program is free software: you can redistribute it and/or modify
% !    it under the terms of the GNU General Public License as published by
% !    the Free Software Foundation, either version 3 of the License, or
% !    (at your option) any later version.
% !
% !    This program is distributed in the hope that it will be useful,
% !    but WITHOUT ANY WARRANTY; without even the implied warranty of
% !    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% !    GNU General Public License for more details.
% !
% !    You should have received a copy of the GNU General Public License
% !    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% !
% !   Contributors list:
% !   - R. Kurnia
% !   - E. van Groesen
% !--------------------------------------------------------------------------------------

fc=fix(clock);

% prev_log=false;
% if exist([Proj.Dir,Proj.savename,'.log'],'file')
% fid = fopen([Proj.Dir,Proj.savename,'.log'],'r');
% textLogPrev = textscan(fid,'%s','delimiter','\n','CollectOutput', true);
% assignin('base','text',textLogPrev)
% fclose(fid);
% delete([Proj.Dir,Proj.savename,'.log']);
% prev_log=true;
% end


diary([Proj.workdir,'Log_',Proj.savename,'.log']);
diary on;
log.string{1}='************************************************************************************';
log.string{2}=['HAWASSI-AB 2D simulation diary, date: ', num2str(fc(3)),'.' num2str(fc(2)),'.' num2str(fc(1)),...
    ' time:  ' num2str(fc(4)),'h.' num2str(fc(5)),'min.'];
log.string{3}='************************************************************************************';
log.string{4}='';
log.string{5}=['Project name: ',Proj.savename];
log.string{6}='';
log.string{7}=['User note   :',Proj.usernote];
log.string{8}='';
log.string{9}='MODEL DESCRIPTION :';
log.string{10}=['Dynamic Model        : ', model.dyn(1:3)];
log.string{11}=['Dispersion Model    : ', model.dispersion(5:end)];
disp(log.string{1})
disp(log.string{2})
disp(log.string{3})
disp('%');
disp(log.string{5});
disp('%');
disp(log.string{7});
disp('%');
disp(log.string{9});
disp(log.string{10});
disp(log.string{11});
idn=11;
if strcmp(model.breaking.check,'Yes')
    log.string{idn+1}=['Breaking             : '];
    log.string{idn+2}=['   Initiation           : U/C= ', num2str(model.breaking.KBC)];
    log.string{idn+3}=['   Termination     : uF/uI= ', num2str(model.breaking.TC),', T*/Tp= ', num2str(model.breaking.Tchar)];
    disp(log.string{idn+1});
    disp(log.string{idn+2});
    disp(log.string{idn+3});
    idn=idn+3;
else
    log.string{idn+1}=['Breaking                    : No'];
    disp(log.string{idn+1});
    idn=idn+1;
end
if model.current.check==1
    log.string{idn+1}=['Current                      :  ux=',num2str(model.current.ux),'[m/s], uy=',num2str(model.current.uy),'[m/s]'];
    disp(log.string{idn+1});
    idn=idn+1;
else
    log.string{idn+1}=['Current                      : No'];
    disp(log.string{idn+1});
    idn=idn+1;
end
log.string{idn+1}='%';
disp(log.string{idn+1});
log.string{idn+2}='INFLUX DESCRIPTION : ';
disp(log.string{idn+2});
idn=idn+2;
if strcmpi(input.wave.option,'No')
    log.string{idn+1}=['There is no influxed wave'];
    disp(log.string{idn+1});
    idn=idn+1;
end

if strcmpi(input.wave.option,'Yes')
    log.string{idn+1}=['Number of influxing        : ', num2str(input.wave.N)];
    disp(log.string{idn+1});
    idn=idn+1;
    for ii=1:input.wave.N
        if input.wave.N>1
        log.string{idn+1}=['--------------------- Influxing-',num2str(ii),'---------------------'];
        disp(log.string{idn+1});
        idn=idn+1;
        end
        if strcmpi(input.wave.type(ii).name,'Harmonic')...
                ||strcmpi(input.wave.type(ii).name,'User-defined (signal)')
            log.string{idn+1}=['Signal type                     : ',cell2mat(input.wave.type(ii).name)];
            disp(log.string{idn+1});
            idn=idn+1;
        elseif strcmpi(input.wave.type(ii).name,'Jonswap')
            log.string{idn+1}=['Signal type                     : ',cell2mat(input.wave.type(ii).name)];
            log.string{idn+2}=['       gamma                      : ',num2str(input.wave.JS_gamma(ii))];
            log.string{idn+3}=['       spreading factor (s): ',num2str(input.wave.JS_s(ii))];
            log.string{idn+4}=['       standard deviation  : ',num2str(rad2deg(influx.Spect(ii).prop.StdTheta),4),' deg'];
            disp(log.string{idn+1});
            disp(log.string{idn+2});
            disp(log.string{idn+3});
            disp(log.string{idn+4});
            idn=idn+4;
        elseif strcmpi(input.wave.type(ii).name,'User-defined (variance density spectrum)')
            log.string{idn+1}=['Signal type                     : ',cell2mat(input.wave.type(ii).name)];
            log.string{idn+2}=['       spreading factor (s): ',num2str(input.wave.userspectrum(ii).spreading)];
            log.string{idn+3}=['       standard deviation  : ',num2str(rad2deg(influx.Spect(ii).prop.StdTheta),4),' deg'];
            disp(log.string{idn+1});
            disp(log.string{idn+2});
            disp(log.string{idn+3});
            idn=idn+3;
        end
        
        if strcmp(input.wave.type(ii).name,'Harmonic')
            log.string{idn+1}=['Amplitude (A)                 : ',num2str(roundn(influx.par.Hs(ii)/2,-3)),'[m]'];
            disp(log.string{idn+1});
            idn=idn+1;
        else
            log.string{idn+1}=['Significant wave Height (Hs)  : ',num2str(roundn(influx.par.Hs(ii),-3)),'[m]'];
            disp(log.string{idn+1});
            idn=idn+1;
        end
        log.string{idn+1}=['Peak period (Tp)           : ',num2str(roundn(influx.par.T_p(ii),-3)),'[s]'];
        disp(log.string{idn+1});
        log.string{idn+2}=['Mean period (Tm01)     : ',num2str(roundn(influx.par.Tm01(ii),-3)),'[s]'];
        disp(log.string{idn+2});
        idn=idn+2;    
    
            if input.wave.type(ii).flag==2 || input.wave.type(ii).flag==3
                    log.string{idn+1}=['Direction       (degree) : ',num2str(roundn(input.wave.direction(ii),-3)),'[deg]'];
                    disp(log.string{idn+1});idn=idn+1; 
            elseif input.wave.type(ii).flag==5
                log.string{idn+1}=['Direction       (degree)  : ',num2str(roundn(input.wave.userspectrum(ii).direction,-3)),'[deg]'];
                disp(log.string{idn+1});idn=idn+1;
            end
        
        tempCat=influx.par.category(ii);
        log.string{idn+1}=['Depth at influxing (h)    : ', num2str(influx.par.meandepth(ii)),'[m]'];
        log.string{idn+2}=['      Derived info:'];
        log.string{idn+3}=['      Peak frequency (nu)                    : ',num2str(roundn(influx.par.nu_p(ii),-3)),'[rad/s]'];
        log.string{idn+4}=['      Mean frequency                           : ',num2str(roundn(2*pi/influx.par.Tm01(ii),-3)),'[rad/s]'];
        log.string{idn+5}=['      Peak wave-number (kp)              : ',num2str(roundn(influx.par.k_p(ii),-3))];
        log.string{idn+6}=['      Peak wave-length                        : ',num2str(roundn(influx.par.lambda_p(ii),-3)),'[m]'];
        log.string{idn+7}=['      Peak phase speed                      : ',num2str(roundn(influx.par.Cp_p(ii),-3)),'[m/s]'];
        log.string{idn+8}=['      Peak group speed                      : ',num2str(roundn(influx.par.Cg_p(ii),-3)),'[m/s]'];
        log.string{idn+9}=['      Steepness (kp*(Hs./2))                : ',num2str(roundn(influx.par.ka(ii),-3))];
        log.string{idn+10}=['      Relative wave-length(lambda/h) : ',num2str(influx.par.lambda_per_D(ii)),'(',tempCat{1:end},')'];
        log.string{idn+11}=['                                                 (kp*h)  : ',num2str(roundn(influx.par.kh(ii),-3))];
        log.string{idn+12}='--------------------------------------------------------------------';
        disp(log.string{idn+1});
        disp(log.string{idn+2});
        disp(log.string{idn+3});
        disp(log.string{idn+4});
        disp(log.string{idn+5});
        disp(log.string{idn+6});
        disp(log.string{idn+7});
        disp(log.string{idn+8});
        disp(log.string{idn+9});
        disp(log.string{idn+10});
        disp(log.string{idn+11});
        disp(log.string{idn+12});
        idn=idn+12;
    end
end
log.string{idn+1}='%';
log.string{idn+2}='INITIAL WAVE CONDITIONS : ';
log.string{idn+3}=['Initial condition        : ',ivp.typename];
disp(log.string{idn+1});
disp(log.string{idn+2});
disp(log.string{idn+3});
idn=idn+3;
if  ivp.type==2||ivp.type==3||ivp.type==4||ivp.type==5
    log.string{idn+1}= ['Amplitude                  : ',num2str(ivp.A)];
    log.string{idn+2}=['Standard deviation : ',num2str(ivp.sigma)];
    log.string{idn+3}=['Center position        : ',num2str(ivp.x0)];
    disp(log.string{idn+1});
    disp(log.string{idn+2});
    disp(log.string{idn+3});idn=idn+3;
end
if ivp.type==6
    log.string{idn+1}=['File name             : ',ivp.userdata.name];
    disp(log.string{idn+1});idn=idn+1;
end

log.string{idn+1}='%';
log.string{idn+2}='BOUNDARY ASSIMILATION : ';
log.string{idn+3}=['Boundary assimilation      : ',input.bdyassim.check];
disp(log.string{idn+1});
disp(log.string{idn+2});
disp(log.string{idn+3});
idn=idn+3;
if input.bdyassim.option==1
    if input.bdyassim.shapeOpt==1, shapeOpt='Half-circle';
    else shapeOpt='User-defined';end;
        
 log.string{idn+1}=['Shape      : ',shapeOpt];
 disp(log.string{idn+1});idn=idn+1;
 if input.bdyassim.shapeOpt==1
 log.string{idn+1}=['at radius:',num2str(input.bdyassim.halfcirc_R1),'; Smooth factor: ',num2str(input.bdyassim.smoothfact)];  
 log.string{idn+2}=['Center position (x,y)=(',num2str(input.bdyassim.halfcirc_xc),...
     ',',num2str(input.bdyassim.halfcirc_yc),') [m]'];  
 disp(log.string{idn+1});
 disp(log.string{idn+2});idn=idn+2;
 end
 log.string{idn+1}=['Propagation direction: ',bdyassim.propdir];
 disp(log.string{idn+1});
 idn=idn+1;  
 tempCat=bdyassim.par.category;
 log.string{idn+1}='Wave info:';
 log.string{idn+2}=['Significant wave Height (Hs)  : ',num2str(roundn(bdyassim.par.Hs,-3)),'[m]'];  
 log.string{idn+3}=['Peak period (Tp)           : ',num2str(roundn(bdyassim.par.Tp,-3)),'[s]'];
 log.string{idn+4}=['      Derived info:'];
 log.string{idn+5}=['      Peak frequency (nu)                    : ',num2str(roundn(bdyassim.par.nupeak,-3)),'[rad/s]'];
 log.string{idn+6}=['      Peak wave-number (kp)              : ',num2str(roundn(bdyassim.par.k_p,-3))];
 log.string{idn+7}=['      Peak wave-length                        : ',num2str(roundn(bdyassim.par.lambda_p,-3)),'[m]'];
 log.string{idn+8}=['      Peak phase speed                      : ',num2str(roundn(bdyassim.par.Cp_p,-3)),'[m/s]'];
 log.string{idn+9}=['      Peak group speed                      : ',num2str(roundn(bdyassim.par.Cg_p,-3)),'[m/s]'];
 log.string{idn+10}=['      Steepness (kp*(Hs./2))                : ',num2str(roundn(bdyassim.par.ka,-3))];
 log.string{idn+11}=['     Relative wave-length(lambda/h) : ',num2str(roundn(bdyassim.par.lambda_p./bdyassim.par.meandepth,-3)),'(',tempCat{1:end},')'];
 log.string{idn+12}=['                                                 (kp*h)  : ',num2str(roundn(bdyassim.par.kh,-3))];
 log.string{idn+13}=['Time interv.  : t start: ', num2str(input.bdyassim.tinit),...
             ' [s]; t end: ', num2str(input.bdyassim.tend),' [s]; time step:',...
             num2str(input.bdyassim.dt),' [s]'];
 log.string{idn+14}='--------------------------------------------------------------------';

 disp(log.string{idn+1});
 disp(log.string{idn+2});
 disp(log.string{idn+3});
 disp(log.string{idn+4});
 disp(log.string{idn+5});
 disp(log.string{idn+6});
 disp(log.string{idn+7});
 disp(log.string{idn+8});
 disp(log.string{idn+9});
 disp(log.string{idn+10});
 disp(log.string{idn+11});
 disp(log.string{idn+12});
 disp(log.string{idn+13});
 disp(log.string{idn+14});
 idn=idn+14;

end

log.string{idn+1}='%';
log.string{idn+2}='NUMERICAL SETTINGS : ';
log.string{idn+3}=['Spatial interval      : ','x in (', num2str(dom.X(1)),';',num2str(dom.X(end)),') [m]'];
log.string{idn+4}=['                                   ','y in (', num2str(dom.Y(1)),';',num2str(dom.Y(end)),') [m]'];
log.string{idn+5}=['Length of Fourier bdy : ','left: ',num2str(dom.fbl.l),'[m]','; right: ',num2str(dom.fbl.r),'[m]'];
log.string{idn+6}=['                                            bottom: ',num2str(dom.fbl.b),'[m]','; top: ',num2str(dom.fbl.t),'[m]'];
log.string{idn+7}=['Number of Nodes : ',' Npx=',num2str(dom.Nx),' Npy=',num2str(dom.Ny),' total:',num2str(dom.Nx*dom.Ny)];
log.string{idn+8}=['        Grid size        : dx= ',num2str(roundn(dom.dx,-3)),'[m]', 'dy= ',num2str(roundn(dom.dy,-3)),'[m]'];
log.string{idn+9}=['        Cutfrac k        : kx: ',num2str(spatial.dom.cutfracx),'; ky: ',num2str(spatial.dom.cutfracy)];
log.string{idn+10}=['Time interval         : ','(', num2str(timeSimul.t_init) ,';',num2str(timeSimul.t_end),') [s]'];
log.string{idn+11}=['        Time step (dt): ', num2str(roundn(timeSimul.dt,-3)),' [s]'];
log.string{idn+12}='%';
disp(log.string{idn+1});
disp(log.string{idn+2});
disp(log.string{idn+3});
disp(log.string{idn+4});
disp(log.string{idn+5});
disp(log.string{idn+6});
disp(log.string{idn+7});
disp(log.string{idn+8});
disp(log.string{idn+9});
disp(log.string{idn+10});
disp(log.string{idn+11});
disp(log.string{idn+12});idn=idn+12;

if strcmpi(bath.type,'Flat')
    log.string{idn+1}=['Bathymetry            : ', bath.name];
    log.string{idn+2}=['     Depth                 : ', num2str(bath.depth),'[m]'];
    disp(log.string{idn+1});
    disp(log.string{idn+2});idn=idn+2;
elseif strcmpi(bath.type,'Slope in (x-axis)') || strcmp(bath.type,'Slope in (y-axis)')
    log.string{idn+1}=['Bathymetry            : ', bath.name];
    log.string{idn+2}=['     Depth            : max= ', num2str(bath.par(1)),'[m],',' min= ', num2str(bath.par(2)),'[m]'];
    if bath.interp==2
    log.string{idn+4}=['Interp method         :', num2str(bath.interp)];
    else
    log.string{idn+4}=['Interp method         :', num2str(bath.interp),'; mid Depth: ', num2str(bath.depthref),'[m]'];
    end
    log.string{idn+3}=['     Slope            : gradient= ', num2str(bath.par(3)), ', start position= ', num2str(bath.par(4)),'[m]'];
    disp(log.string{idn+1})
    disp(log.string{idn+2});
    disp(log.string{idn+3});
    disp(log.string{idn+4});idn=idn+4;
elseif strcmp(bath.type,'Shore (in x-axis)') || strcmp(bath.type,'Shore (in y-axis)')
    log.string{idn+1}=['Bathymetry            : ', bath.name];
    log.string{idn+2}=['     Depth            : max= ', num2str(bath.par(1)),'[m]; min total depth=',num2str(bath.Hmin),'[m]'];
    log.string{idn+3}=['     Slope            : gradient=', num2str(bath.par(2)), ', shore position : ', num2str(bath.par(3)),'[m]'];
    log.string{idn+4}=['Interp method         :', num2str(bath.interp),'; mid Depth: ', num2str(bath.depthref),'[m]'];
    disp(log.string{idn+1})
    disp(log.string{idn+2});
    disp(log.string{idn+3});
    disp(log.string{idn+4});idn=idn+4;
else
    log.string{idn+1}=['Bathymetry            : ', bath.name];
    if bath.interp==2
    log.string{idn+2}=['Interp method         :', num2str(bath.interp)];
    else
    log.string{idn+2}=['Interp method         :', num2str(bath.interp),'; mid Depth: ', num2str(bath.depthref),'[m]'];
    end
    disp(log.string{idn+1});
    disp(log.string{idn+2});idn=idn+2;
    if strcmpi(bath.name,'Shore')
    log.string{idn+1}=['A minimum total depth for shoreline: ',num2str(bath.Hmin),'[m]'];
    disp(log.string{idn+1});idn=idn+1;
    end
end

if bath.friction.check==1
    log.string{idn+1}='%';
    disp(log.string{idn+1});idn=idn+1;
    fricdata=bath.friction.param;
    Ni=length(fricdata(:,1));
    log.string{idn+1}='Bottom friction:';
    disp(log.string{idn+1});idn=idn+1;
    for ii=1:Ni
    if strcmpi(fricdata(ii,1),'Rectangle')
    log.string{idn+1}=['#', num2str(ii), ' Position: x in (',num2str(cell2mat(fricdata(ii,2))),';',...
        num2str(cell2mat(fricdata(ii,3))),'),', ' y in (',num2str(cell2mat(fricdata(ii,4))),';',...
        num2str(cell2mat(fricdata(ii,5))),');', ' Friction Coef.= ',num2str(cell2mat(fricdata(ii,6)))];
    else
    log.string{idn+1}=['#', num2str(ii), ' User-defined;', ' Friction Coef.= ',num2str(cell2mat(fricdata(ii,6)))];    
    end
    disp(log.string{idn+1});idn=idn+1;
    end
end

if strcmpi(input.wall.option,'Yes')
    log.string{idn+1}='%';
    disp(log.string{idn+1});idn=idn+1;
    wallparam=dom.wall.param;
    log.string{idn+1}=['Number of wall: ', num2str(length(wallparam(:,1)))];
    disp(log.string{idn+1});idn=idn+1;
    for ii=1:length(wallparam(:,1))
    log.string{idn+1}=['#',num2str(ii),' shape  :', cell2mat(wallparam(ii,1)),'; Method  :', cell2mat(wallparam(ii,9)),'; Type  :', cell2mat(wallparam(ii,10)),'; Refl. coef   :', num2str(cell2mat(wallparam(ii,11)))];
    disp(log.string{idn+1});idn=idn+1;   
    end
end
log.string{idn+1}='%';
disp(log.string{idn+1});idn=idn+1;


if strcmpi(input.wave.option,'Yes')
    log.string{idn+1}=['Influxing setting'];
    disp(log.string{idn+1});idn=idn+1;
    
    for ii=1:input.wave.N
        if input.wave.N>1
        log.string{idn+1}=['--------------------- Influxing-',num2str(ii),'--------------------------------'];
        disp(log.string{idn+1});idn=idn+1;    
        end
        log.string{idn+1}=['Influx method : ', cell2mat(model.influx.type(ii)),' (', cell2mat(model.influx.direction(ii)),')'];
        disp(log.string{idn+1});idn=idn+1;
        if strcmpi(spatial.influx.linetype(ii),'Straight')
          log.string{idn+1}=['Straight line : start point at (x,y)=(', num2str(cell2mat(spatial.influx.x1(ii))),...
             ',',num2str(cell2mat(spatial.influx.y1(ii))),')[m]; end point at (x,y)=(', num2str(cell2mat(spatial.influx.x2(ii))),',',...
             num2str(cell2mat(spatial.influx.y2(ii))),') [m]'];
         disp(log.string{idn+1});idn=idn+1;
        else
          log.string{idn+1}=['Arc circle line : center point at (x,y)= (', num2str(cell2mat(spatial.influx.xc(ii))),...
             ',',num2str(cell2mat(spatial.influx.yc(ii))),') [m], radius=', num2str(cell2mat(spatial.influx.radius(ii))),'[m], with angle in (',...
             num2str(cell2mat(spatial.influx.theta1(ii))),',',num2str(cell2mat(spatial.influx.theta2(ii))),') deg'];
          disp(log.string{idn+1});idn=idn+1;    
        end
        log.string{idn+1}=['Time interv.  : t start: ', num2str(input.wave.t_init(ii)),...
             ' [s]; t end: ', num2str(input.wave.t_end(ii)),' [s]; time step:',...
             num2str(input.wave.dt(ii)),' [s]'];
        log.string{idn+2}='--------------------------------------------------------------------';
         disp(log.string{idn+1}) 
         disp(log.string{idn+2});idn=idn+2;     
    end
end
log.string{idn+1}='%';
disp(log.string{idn+1});idn=idn+1;

% if ~strcmp(model.influx.direction,'None')
%     if strcmp(model.influx.type,'Point')
%         disp(['Influx method : ', model.influx.direction{1:end},'; ', model.influx.type{1:end}]);
%     else
%         disp(['Influx method : ', model.influx.direction{1:end},'; ', model.influx.type{1:end}]);
%     end
%     
%     if input.wave.adjzone.check==1
%         disp(['    Nonlinear adjustment : ', num2str(input.wave.adjzone.lengthfact), '*lambda_peak']);
%     else
%         disp(['    Nonlinear adjustment : 0*lambda_peak']);
%     end
%     if input.wave.ramp.check==1
%         disp(['    Ramp                 : ', num2str(input.wave.ramp.length), '*Tp']);
%     else
%         disp(['    Ramp                 : Off']);
%     end
%     if input.wave.filter.check==1
%         disp(['    Filter               : On; Freq. interval: (',num2str(input.wave.filter.LFreq),';',num2str(input.wave.filter.HFreq),') rad/s']);
%     else
%         disp(['    Filter               : Off']);
%     end
%     
%     disp(['    Influx position      : x in (', num2str(spatial.influx.x(1)),',',num2str(spatial.influx.x(2)), ') [m]']);
%     disp(['                           y in (', num2str(spatial.influx.y(1)),',',num2str(spatial.influx.y(2)), ') [m]']);
% else
%     disp(['Influx method : ', model.influx.direction{1:end}]);
% end
% disp('%');


% if options.interior.check==1
%     disp(['INTERIOR CALCULATION PREPARATION: ']);
%     disp(['Time interval        : (',num2str(options.interior.time(1)),';',num2str(options.interior.time(2)),') [s]']);
%     disp(['Time step            : ',num2str(par.dt.*options.interior.time(3)),' [s]']);
%     disp('%');
% end
log.string{idn+1}='ODE SETTING:';
log.string{idn+2}=['ODE solver           : ', par.ode.solver];
log.string{idn+3}=['Relative error tol.  : ', num2str(par.ode.tol)];
log.string{idn+4}=['No of Partition      : ', num2str(par.ode.Npartition)];
log.string{idn+5}='---------------------------------------------------------------------';
disp(log.string{idn+1});
disp(log.string{idn+2});
disp(log.string{idn+3});
disp(log.string{idn+4});
disp(log.string{idn+5});
diary off
%%diary off; in HaWaSSI.m in function of numerical_setting_advises

%     if prev_log
%     disp(textLogPrev{1,1});
%     end
% if ~isdeployed
% clc;
% end
