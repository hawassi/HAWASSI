function plot_PP1(CheckId,X,Z,P_t,Pdyn_t,Pnonlindyn_t,V_t,Vx_t,Vz_t,ax_t,az_t,StringVar,stringlabel,stringtitle,SetAx,Id)
Check_Vz=CheckId.Vz;
Check_Vx=CheckId.Vx;
Check_V=CheckId.V;
Check_ax=CheckId.ax;
Check_az=CheckId.az;
Check_P=CheckId.P;
Check_Pdyn=CheckId.Pdyn;
Check_Pnonlindyn=CheckId.Pnonlindyn;




if Check_Vz==1 && Check_Vx==0 && Check_V==0 && Check_P==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_Vx==1 && Check_V==0 && Check_P==0 && Check_Vz==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_V==1 && Check_P==0 && Check_Vz==0 && Check_Vx==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
%   plot_PP1_V(X,Z,V_t,Vx_t,Vz_t,StringVar{2},stringlabel,SetAx,Id);
    title(stringtitle);
elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    title(stringtitle);
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    
    
elseif Check_P==1  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{5},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1&& Check_Vx==0 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==1 && Check_V==0&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==1&&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
 
elseif Check_P==0  && Check_Vz==0&& Check_Vx==0 && Check_V==0&&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(2,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)

elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)   
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
        title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)   
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)

   elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
     title(stringtitle);   
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)   
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)

    
  
 elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
  elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
  elseif Check_P==0  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
       subplot(3,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
 
    
    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
  elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)   
  elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
     title(stringtitle);
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)   
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
  
   elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(3,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)
    title(stringtitle);   
    subplot(3,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)   
    subplot(3,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
       
    
    
    
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)   
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)   
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)   
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)   
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)    

elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)   
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)   
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)    
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)     
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)  
    
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)     
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)  
    
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)     
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)  
 
elseif Check_P==1  && Check_Vz==0 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)     
    subplot(4,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(4,1,4)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    
elseif Check_P==0  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)    
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
    
elseif Check_P==0  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(4,1,1)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(4,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id)    
    subplot(4,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(4,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)    
    subplot(5,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(5,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==0  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
        plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id)  
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{5},stringlabel,SetAx,Id)    
    subplot(5,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(5,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==0 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
    subplot(5,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(5,1,4)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
elseif Check_P==1  && Check_Vz==1 && Check_Vx==0 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
    subplot(5,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(5,1,4)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)    
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==0 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
    subplot(5,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(5,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id)
  elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==0&&Check_Pnonlindyn==1
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(5,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id) 
    subplot(5,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id) 
    elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==0
    subplot(5,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(5,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
    subplot(5,1,3)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id) 
    subplot(5,1,4)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(5,1,5)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id) 
 elseif Check_P==1  && Check_Vz==1 && Check_Vx==1 && Check_V==1 &&Check_Pdyn==1&&Check_Pnonlindyn==1
    subplot(6,1,1)
    plot_PP1_Var(X,Z,P_t,StringVar{1},stringlabel,SetAx,Id)      
    title(stringtitle);
    subplot(6,1,2)
    plot_PP1_Var(X,Z,Pdyn_t,StringVar{5},stringlabel,SetAx,Id) 
    subplot(6,1,3)
    plot_PP1_Var(X,Z,Pnonlindyn_t,StringVar{6},stringlabel,SetAx,Id) 
    subplot(6,1,4)
    plot_PP1_Var(X,Z,V_t,StringVar{2},stringlabel,SetAx,Id)
    subplot(6,1,5)
    plot_PP1_Var(X,Z,Vx_t,StringVar{3},stringlabel,SetAx,Id)
    subplot(6,1,6)
    plot_PP1_Var(X,Z,Vz_t,StringVar{4},stringlabel,SetAx,Id) 
  else
    
    return;
end

if  Check_ax==1 && Check_az==1
    figure
    subplot(2,1,1)
    plot_PP1_Var(X,Z,ax_t,'dxdtPHI\_CFD',stringlabel,SetAx,Id)
    title(stringtitle);
    subplot(2,1,2)
    plot_PP1_Var(X,Z,az_t,'dtdxPHI\_CFD',stringlabel,SetAx,Id)
end
