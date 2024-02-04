%%%%%%%%%%%%%%%%Internal Flow module Preparation%%%%%%%%%%%%%%%%%%%%%%%%%%%
global   dteta dtphihat iterInterior timeIP uIP etaIP
dteta=[]; dtphihat=[]; iterInterior=1; timeIP=[]; uIP=[]; etaIP=[];

if options.interior.check==1
    
    
    if tintval(1)> par.IP.time(2) || tintval(end)< par.IP.time(1)
        checkInterior=0;
    else
        checkInterior=1;
        if Npartition==1
            t_IP_start=par.IP.time(1);
            t_IP_end  =par.IP.time(2);
        else
            if tintval(1)> par.IP.time(1)
                t_IP_start=tintval(1);
            else
                t_IP_start=par.IP.time(1);
            end
            
            if tintval(end)< par.IP.time(2)
                t_IP_end=tintval(end);
            else
                t_IP_end=par.IP.time(2);
            end
            
        end
        
        
        Nt_IP=round((t_IP_end-t_IP_start)/(par.IP.time(3).*par.dt));
        dteta=zeros(Nt_IP,length(par.x));
        dtphihat=zeros(Nt_IP,length(par.x));
        timeIP  =zeros(Nt_IP,1);
        uIP=zeros(Nt_IP,length(par.x));
        etaIP=zeros(Nt_IP,length(par.x));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
