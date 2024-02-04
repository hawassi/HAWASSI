%%%Input data for post-processing and Internal flows module
global Idstop
        
inputPProc.time=output.time; inputPProc.eta=output.eta;
inputPProc.u=output.u;
if shippar.check== 1
inputPProc.shipRB=output.RB;
inputPProc.shipsavevar=shipsavevar;
end
inputPProc.x=x; inputPProc.bathy=par.bathy;
inputPProc.Wall=par.wall;
inputPProc.dynmodel=model.dyn;
inputPProc.path=Proj.path;
inputPProc.savename=Proj.savename;
inputPProc.Xinflux=par.Xinflux;
inputPProc.dataBreak_nodes=output.break_nodes;
inputPProc.break_crest    =output.break_crest;
inputPProc.Fbdy   =par.bf0;
inputPProc.lambda_p   =influx.lambda_p;
inputPProc.Nonlin_Adj =bath.influx_AdjZone;
if input.type ==3
    inputPProc.tcoarse=influx.tcoarse;
end


if IndtIntF>0
        if  Idstop==1 &&checkInterior==1
         temp_tIP=timeIP(1:IndtIntF,1);
         timeIP=[];
         timeIP=temp_tIP;
         
         temp_dteta=dteta(1:IndtIntF,:);
         dteta=[];
         dteta=temp_dteta;
         
         temp_dtphihat=dtphihat(1:IndtIntF,:);
         dtphihat=[];
         dtphihat=temp_dtphihat;
         
         temp_uIPAll=uIPAll(1:IndtIntF,:);
         uIPAll=[];
         uIPAll=temp_uIPAll;
         
         temp_etaIPAll=etaIPAll(1:IndtIntF,:);
         etaIPAll=[];
         etaIPAll=temp_etaIPAll;
        end
        
        
        InputInteriorCalc.dteta   =dteta;
        InputInteriorCalc.dtphihat=dtphihat;
        Om         = str2func(model.dispersion);
        InputInteriorCalc.Om      =Om;
        InputInteriorCalc.k       =k;
        InputInteriorCalc.nupeak  =influx.nu_p;
        InputInteriorCalc.x       =x;
        InputInteriorCalc.u       =uIPAll;
        InputInteriorCalc.eta     =etaIPAll;
        InputInteriorCalc.time    =timeIP;
        InputInteriorCalc.timeInterv=par.IP.time;
        InputInteriorCalc.bathy    =par.bathy;
        InputInteriorCalc.Fbdy     =par.bf0;
        if any(par.bathy>0) %%Run-up case
            InputInteriorCalc.H_min=H_min;
        end
        InputInteriorCalc.path=Proj.path;
        InputInteriorCalc.savename=Proj.savename;
end

