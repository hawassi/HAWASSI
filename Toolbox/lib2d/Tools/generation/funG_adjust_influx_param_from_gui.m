function [inputwave,modelinflux,spatialinflux]=funG_adjust_influx_param_from_gui(inputwave,dom)
propdata=inputwave.propertiesdata;
methoddata=inputwave.methoddata;
userdata=inputwave.userdata;

modelinflux.type=methoddata(:,1);
spatialinflux.linetype  =methoddata(:,2);
spatialinflux.x1    =methoddata(:,3);
spatialinflux.y1    =methoddata(:,4);
spatialinflux.x2    =methoddata(:,5);
spatialinflux.y2    =methoddata(:,6);
spatialinflux.xc    =methoddata(:,7);
spatialinflux.yc    =methoddata(:,8);
spatialinflux.radius=methoddata(:,9);
spatialinflux.theta1=methoddata(:,10);
spatialinflux.theta2=methoddata(:,11);

inputwave.N=length(propdata(:,1));
for ii=1:inputwave.N
    spatialinflux.line(ii).Orientation=funG_lineOrientation(spatialinflux,ii);
    spatialinflux.line(ii).xy=funG_lineXYmapped(dom,spatialinflux,ii);
    if strcmpi(propdata(ii,1),'Harmonic')
        inputwave.type(ii).flag=2;
        if iscellstr(propdata(ii,2))
            propdata(ii,2)={str2num(cell2mat(propdata(ii,2)))};
        end
        if iscellstr(propdata(ii,4))
            propdata(ii,4)={str2num(cell2mat(propdata(ii,4)))};
        end
        if iscellstr(propdata(ii,7))
            propdata(ii,7)={str2num(cell2mat(propdata(ii,7)))};
        end
        inputwave.Hs(ii)=cell2mat(propdata(ii,2))*2;
        inputwave.Tp(ii)=cell2mat(propdata(ii,4));
        inputwave.direction(ii)=cell2mat(propdata(ii,7));
        
        modelinflux.direction(ii)=funG_findInfluxDir(spatialinflux.line(ii).Orientation,inputwave.direction(ii));
        inputwave.type(ii).name={'Harmonic'};
    elseif strcmpi(propdata(ii,1),'JONSWAP')
        if iscellstr(propdata(ii,3))
            propdata(ii,3)={str2num(cell2mat(propdata(ii,3)))};
        end
        if iscellstr(propdata(ii,4))
            propdata(ii,4)={str2num(cell2mat(propdata(ii,4)))};
        end
        if iscellstr(propdata(ii,5))
            propdata(ii,5)={str2num(cell2mat(propdata(ii,5)))};
        end
        if iscellstr(propdata(ii,6))
            propdata(ii,6)={str2num(cell2mat(propdata(ii,6)))};
        end
        if iscellstr(propdata(ii,7))
            propdata(ii,7)={str2num(cell2mat(propdata(ii,7)))};
        end
        
        inputwave.type(ii).flag=3;
        inputwave.Hs(ii)=cell2mat(propdata(ii,3));
        inputwave.Tp(ii)=cell2mat(propdata(ii,4));
        inputwave.JS_gamma(ii)=cell2mat(propdata(ii,5));
        inputwave.JS_s(ii)=cell2mat(propdata(ii,6));
        inputwave.direction(ii)=cell2mat(propdata(ii,7));
        
        modelinflux.direction(ii)=funG_findInfluxDir(spatialinflux.line(ii).Orientation,inputwave.direction(ii));
        inputwave.type(ii).name={'JONSWAP'};
    elseif strcmpi(propdata(ii,1),'User-defined (signal)')
        inputwave.type(ii).flag=4;
        inputwave.usersignal(ii).inflX=userdata(ii).inflX;
        inputwave.usersignal(ii).inflY=userdata(ii).inflY;
        inputwave.usersignal(ii).time=userdata(ii).time;
        inputwave.usersignal(ii).eta=userdata(ii).eta;
        modelinflux.direction(ii)=funG_findInfluxDirUserDefined(spatialinflux,ii,dom);
        inputwave.type(ii).name={'User-defined (signal)'};
    elseif strcmpi(propdata(ii,1),'User-defined (variance density spectrum)')
        if iscellstr(propdata(ii,6))
            propdata(ii,6)={str2num(cell2mat(propdata(ii,6)))};
        end
        if iscellstr(propdata(ii,7))
            inputwave.direction(ii)={str2num(cell2mat(propdata(ii,7)))};
        end
        
        inputwave.type(ii).flag=5;
        inputwave.userspectrum(ii).varDensity=userdata(ii).varianceDensity(:,2);
        inputwave.userspectrum(ii).ww=userdata(ii).varianceDensity(:,1);
        inputwave.userspectrum(ii).spreading=cell2mat(propdata(ii,6));
        inputwave.userspectrum(ii).direction=cell2mat(propdata(ii,7));
        modelinflux.direction(ii)=funG_findInfluxDir(spatialinflux.line(ii).Orientation(ii),inputwave.userspectrum(ii).direction);
        inputwave.type(ii).name={'User-defined (variance density spectrum)'};
    end
    
    if iscellstr(methoddata(ii,12))
        methoddata(ii,11)={str2num(cell2mat(methoddata(ii,12)))};
    end
    if iscellstr(methoddata(ii,13))
        methoddata(ii,12)={str2num(cell2mat(methoddata(ii,13)))};
    end
    if iscellstr(methoddata(ii,14))
        methoddata(ii,13)={str2num(cell2mat(methoddata(ii,14)))};
    end
    inputwave.t_init(ii)=cell2mat(methoddata(ii,12));
    inputwave.t_end(ii)=cell2mat(methoddata(ii,13));
    inputwave.dt(ii)=cell2mat(methoddata(ii,14));
end