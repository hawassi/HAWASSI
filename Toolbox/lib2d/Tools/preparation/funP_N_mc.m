function N=funP_N_mc(mc)
if mc.check==1
    N=1;
    data=mc.waveinput_data;
    if mc.waveinput_check==1
      for i=1:length(data)  
          if isfield(data(i),'A_param')
          N=max(N,length(data(i).A_param));
          end
          if isfield(data(i),'Tp_param')
          N=max(N,length(data(i).Tp_param));
          end
          if isfield(data(i),'Hs_param')
          N=max(N,length(data(i).Hs_param));
          end
          if isfield(data(i),'gamma_param')
          N=max(N,length(data(i).gamma_param));
          end
          if isfield(data(i),'s_param')
          N=max(N,length(data(i).s_param));
          end
      end
    end
    
    if mc.numset_check==1
       N=max(N,length(mc.numset_px));
       N=max(N,length(mc.numset_py));
    end
    
    if mc.numrun_check==1
       N=mc.numrun-1;
    end
    
else
   N=1; 
end