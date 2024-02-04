%%%Input data for post-processing and Internal flows module

if input.assim.check==1
 inputPProc.time=Time; inputPProc.eta=Eta;
inputPProc.u=U;   
else   
inputPProc.time=time; inputPProc.eta=eta;
inputPProc.u=u;
end

if isfield(par,'interp')
 inputPProc.OprtInterp_par=par.interp; 
end
inputPProc.x=x; inputPProc.bathy=par.bathy;
inputPProc.Wall=par.wall;
inputPProc.dynmodel=model.dyn;
inputPProc.Proj=Proj;
inputPProc.savename=Proj.savename;
inputPProc.Xinflux=par.Xinflux;
inputPProc.influx=influx;
inputPProc.dataBreak_nodes=dataBreak_nodes;
inputPProc.break_crest    =dataCrestBreak;
inputPProc.Fbdy   =par.bf0;
inputPProc.lambda_p   =influx.lambda_p;
inputPProc.Nonlin_Adj =bath.influx_AdjZone;
inputPProc.bath    =bath;
inputPProc.Oprt    =Oprt;
inputPProc.shippar=shippar;
if shippar.check==1
inputPProc.shipRB=RB; 
inputPProc.shipsavevar=shipsavevar;
end

if input.type ==4
    inputPProc.tcoarse=influx.tcoarse;
end
    
     
global Idstop   
if checkInterior==1        
    if time(end)<par.IP.time(1)
        checkInterior=0;
    else
       if  Idstop==1   &&time(end)< par.IP.time(2)  
         ind0=find(timeIP(2:end)==0,1,'first')+1;
         temp_t=timeIP(1:ind0-1,1);
         timeIP=[];
         timeIP=temp_t;
         
         temp_dteta=dteta(1:ind0-1,:);
         dteta=[];
         dteta=temp_dteta;
         
         temp_dtphihat=dtphihat(1:ind0-1,:);
         dtphihat=[];
         dtphihat=temp_dtphihat;
         
         temp_uIP=uIP(1:ind0-1,:);
         uIP=[];
         uIP=temp_uIP;
         
         temp_etaIP=etaIP(1:ind0-1,:);
         etaIP=[];
         etaIP=temp_etaIP;
        end
       

        InputInteriorCalc.dteta   =dteta;
        InputInteriorCalc.dtphihat=dtphihat;
        
        Om         = str2func(model.dispersion);
        InputInteriorCalc.Om      =Om;
        InputInteriorCalc.k       =k;
        InputInteriorCalc.nupeak  =influx.nu_p;
        InputInteriorCalc.x       =x;
        InputInteriorCalc.u       =uIP;
        InputInteriorCalc.eta     =etaIP;
        InputInteriorCalc.time    =timeIP;
        InputInteriorCalc.timeInterv=par.IP.time;
        InputInteriorCalc.bathy    =par.bathy;
        InputInteriorCalc.Fbdy     =par.bf0;
        if any(par.bathy>0) %%Run-up case
            InputInteriorCalc.H_min=par.interp.H_min;
        end
        InputInteriorCalc.Proj=Proj;
        InputInteriorCalc.savename=Proj.savename;
        if Npartition==1
            save ('-v7.3',[Proj.Dir,Proj.savename,'_simul_Interior.mat'],'InputInteriorCalc')
        else
            inputPProc.simul_file=[Proj.savename,'_simul_part_',num2str(I),'.mat'];
            save ('-v7.3',[Proj.Dir,Proj.savename,'_simul_Interior_part_',num2str(I),'.mat'],'InputInteriorCalc')
        end
    end
end

