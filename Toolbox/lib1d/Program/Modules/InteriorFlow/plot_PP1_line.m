function plot_PP1_line(CheckId,X,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,StringVar,stringlabel,stringtitle)
Check_Vz=CheckId.Vz;
Check_Vx=CheckId.Vx;
Check_V=CheckId.V;
Check_P=CheckId.P;
Check_Pdyn=CheckId.Pdyn;
Check_Pnonlindyn=CheckId.Pnonlindyn;

if  Check_Vz==1 && Check_Vx==0 && Check_V==0 && Check_P==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    title(stringtitle);
elseif Check_Vx==1 && Check_V==0 && Check_P==0 && Check_Vz==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    title(stringtitle);
elseif Check_V==1 && Check_P==0 && Check_Vz==0 && Check_Vx==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    title(stringtitle);
elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    plot_PP1_Var_line(X,Pdyn_t,StringVar{1},stringlabel)
    title(stringtitle);
elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{1},stringlabel)
    title(stringtitle);
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{2},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    
    
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{2},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{5},stringlabel)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
 
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)

elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)  
    subplot(3,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    title(stringtitle);
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)   
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)   
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    title(stringtitle);
   elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)   
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    title(stringtitle);   
    
  
 elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
  elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
  elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
       subplot(3,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
 
    
    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(3,1,2)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    title(stringtitle);
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)   
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)   
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    title(stringtitle);
   elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)
    subplot(3,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)   
    subplot(3,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    title(stringtitle);      
    
    
    
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,3)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)   
    subplot(4,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)   
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)   
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)   
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)    

elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)   
    subplot(4,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,3)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)   
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)    
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)     
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)  
    
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)     
    subplot(4,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)  
    
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)     
    subplot(4,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)  
 
elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)     
    subplot(4,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(4,1,4)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    
elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)    
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
    
elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel)    
    subplot(4,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(4,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)  
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)    
    subplot(5,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(5,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
        plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel)  
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{5},stringlabel)    
    subplot(5,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(5,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(5,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(5,1,4)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(5,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(5,1,4)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(5,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(5,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel)
  elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(5,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel) 
    subplot(5,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel) 
    elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(5,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(5,1,3)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel) 
    subplot(5,1,4)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(5,1,5)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel) 
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(6,1,1)
    plot_PP1_Var_line(X,P_t,StringVar{1},stringlabel)      
    title(stringtitle);
    subplot(6,1,2)
    plot_PP1_Var_line(X,Pdyn_t,StringVar{5},stringlabel) 
    subplot(6,1,3)
    plot_PP1_Var_line(X,Pnonlindyn_t,StringVar{6},stringlabel) 
    subplot(6,1,4)
    plot_PP1_Var_line(X,V_t,StringVar{2},stringlabel)
    subplot(6,1,5)
    plot_PP1_Var_line(X,Vx_t,StringVar{3},stringlabel)
    subplot(6,1,6)
    plot_PP1_Var_line(X,Vz_t,StringVar{4},stringlabel) 
  else
    
    return;
end