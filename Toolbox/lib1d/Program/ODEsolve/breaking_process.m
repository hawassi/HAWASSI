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

%Version 25 Oct 2014
function [B] = breaking_process(eta,u,depth,t,x,param,indEND,t_init)
n       = length(eta);
U       = abs(u);
T_star  = param.Tstar;%[Tp/10;Tp/2]
b       = param.KBC;%Upart>=b*Ccrest, N_I=b*Ccrest
d       = param.TC;%N_F=d* N_I;
B       = zeros(n,1);

global flag iterNprev  CrestBreakPrev CrestBreakNew ...
    add indxMaxNow ID_Node_end  NCBprev ITERbn...
    dataBreak_nodes ITERbdt Break_nodes Node_end  dataCrestBreak

%% Kinematic criteria U>C and find all peaks
[CrestBreak] = Kinematic_breaking_criterion(t,eta,u,x,depth,b,param,indEND);

if ~isempty(CrestBreak)&&flag==0
    flag=1;
    CrestBreakNew=CrestBreak;
elseif isempty(CrestBreak) &&flag==2
    CrestBreakNew=CrestBreakPrev;
    NCBprev=length(CrestBreakPrev.index);
elseif ~isempty(CrestBreak)&& flag==2
    %%
    NCB=length(CrestBreak.index) ;
    NCBprev=length(CrestBreakPrev.index);
    CrestBreakNew=CrestBreakPrev;
    
    add=1;
    for k=1:NCB
        IdCheckNewOld=1;
        for l=1:NCBprev
            dleft=abs(ID_Node_end(l)-CrestBreakPrev.indexNow(l));            
            if  CrestBreakPrev.DirProp(l)==1
                if (CrestBreak.index(k)>=CrestBreakPrev.indexNow(l)-3*dleft) && (CrestBreak.index(k)<=ID_Node_end(l))
                    IdCheckNewOld=0; %restrict to different wave
                    break;
                end
            else
                if (CrestBreak.index(k)>=CrestBreakPrev.indexNow(l)) && (CrestBreak.index(k)<=ID_Node_end(l)+3*dleft)
                    IdCheckNewOld=0; %restrict to different wave
                    break;
                end                
            end
        end
        if IdCheckNewOld ==1 % It will be saved if Crest break candidates are not in the previous breaking nodes 
            CrestBreakNew.index(NCBprev+add)   =CrestBreak.index(k);%|
            CrestBreakNew.indexNow(NCBprev+add)=CrestBreak.index(k);%|same for the first time
            ID_Node_end(NCBprev+add)           =CrestBreak.index(k);%| 
            CrestBreakNew.Ucrest(NCBprev+add)  =CrestBreak.Ucrest(k);
            CrestBreakNew.position(NCBprev+add)=CrestBreak.position(k);
            CrestBreakNew.time(NCBprev+add)    =CrestBreak.time(k);
            CrestBreakNew.DirProp(NCBprev+add) =CrestBreak.DirProp(k);
            add=add+1;
        end
    end
    
%   NCBnow=length(CrestBreakNew.index);
%   %  if NCBnow>NCBprev
%         disp(['=======Multiple Breaking [',num2str(ITERbdt),']=======']);
%         disp(['Number of Breaking: ',num2str(NCBnow)]);
%         disp(['Crest Positions (Start): ',num2str(x([CrestBreakNew.index]))]);
%         disp(['Crest Positions (Now)  : ',num2str(x([CrestBreakNew.indexNow]))]);
%   %  end
   
end

if (flag==1 || flag==2)
   NCBnow=length(CrestBreakNew.index);% dx=x(2)-x(1);
%     if flag==1
%         if NCBnow>1
%             disp(['=======Breaking[',num2str(ITERbdt),']=======']);
%             disp(['Number of Breaking: ',num2str(NCBnow)]);
%             disp(['Positions: ',num2str(CrestBreakNew.position)]);
%         else
%             disp(['=======Breaking[',num2str(ITERbdt),']======='])
%              disp(['Number of Breaking: ',num2str(NCBnow)]);
%             disp(['Crest Positions (Start): ',num2str(x([CrestBreakNew.index]))]);
%             if isfield(CrestBreakNew,'indexNow')
%             disp(['Crest Positions (Now)  : ',num2str(x([CrestBreakNew.indexNow]))]);
%             end
%         end
%     end
%     
    CrestBreakPrev=[];
    iterNprev=1;
    IDTerminate=zeros(NCBnow,1);
     
    it_indexb=1;Break_nodes=[];
    for lx=1:NCBnow
        DirProp=CrestBreakNew.DirProp(lx);
        
        if flag==1
            indxMaxNow=CrestBreakNew.index(lx);
        else
            indxMaxStart=CrestBreakNew.indexNow(lx);
            indxMaxLoc=closest(eta(indxMaxStart:DirProp:ID_Node_end(lx)),...
                             max(eta(indxMaxStart:DirProp:ID_Node_end(lx))));
            indxMaxNow=indxMaxStart+ (indxMaxLoc-1).*DirProp;
        end
        CrestBreakNew.indexNow(lx)=indxMaxNow;
          
        if lx~=1
         if indxMaxNow== CrestBreakNew.indexNow(lx-1);
         IDTerminate(lx)=0;    % terminate a CrestBreak which is same as previous one;
         continue;
         end
        end
        
        tbreak=CrestBreakNew.time(lx);
        U_I=CrestBreakNew.Ucrest(lx);
        U_F=d.*U_I;
        if  t-tbreak>=T_star
            U_star=U_F;
        elseif (t-tbreak >=0) && (t-tbreak <T_star)
            U_star=U_I+((t-tbreak)/T_star).*(U_F-U_I);
        else
            U_star=U_I;
        end
        
        if isempty(indEND)
            if DirProp==1
                indEND=param.xbend;
            else
                indEND=param.xbstart;
            end
        end

        
        Node_end=indxMaxNow;
        for ly=indxMaxNow:DirProp:indEND %index node from the crest
            U_FF=U_F;
            if   U(ly)< U_FF 
                Node_end=ly;
                break;
            end

            break_index=ly;

            if U(break_index)> 2*U_star
                B(break_index)  =  1;
            elseif U(break_index)<= U_star
                B(break_index)  =  0;
            else
                B(break_index)  = (U(break_index)/ U_star )-1;
            end
            Break_nodes(it_indexb)=break_index;
            it_indexb=it_indexb+1;
        end
        
        
%         xstartxend=[x(indxMaxNow) x(Node_end)]
       
        if all(U(indxMaxNow:DirProp:Node_end)<U_star) %&& t-tbreak>0
            IDTerminate(lx)=0;
%             IndexStart=CrestBreakNew.index(lx);
%             xb=x(IndexStart);
%             tb=CrestBreakNew.time(lx);
% %             display (['xb=',num2str(xb),'[m];',' tb=',num2str(tb),'[s];',' dtb= ',num2str(t-tbreak),'[s].']);
%             databreak(ITER,1)=xb;
%             databreak(ITER,2)=tb;
%             databreak(ITER,3)=t-tb;
%             ITER=ITER+1;
        else
            IDTerminate(lx)=1;
            ID_Node_end(iterNprev)=Node_end;
            CrestBreakPrev.index(iterNprev)=CrestBreakNew.index(lx);
            CrestBreakPrev.indexNow(iterNprev)=CrestBreakNew.indexNow(lx);
            CrestBreakPrev.Ucrest(iterNprev)=CrestBreakNew.Ucrest(lx);
            CrestBreakPrev.position(iterNprev)=CrestBreakNew.position(lx);
            CrestBreakPrev.time(iterNprev)=CrestBreakNew.time(lx);
            CrestBreakPrev.DirProp(iterNprev)=DirProp;
            iterNprev=iterNprev+1;
        end
    end
   
    
    if all(IDTerminate==0)
        flag=0;
    else
        flag=2;
    end
    
    if t>t_init+ITERbdt*param.dt  
%      if ITERbn==1, figure; end
%          subplot(4,1,1)
%          plot(x,eta,'r',x([Break_nodes]),eta([Break_nodes]),'ob',x([Break_nodes]),B([Break_nodes]),'*g')
%          title(['time= ', num2str(t)])
%          ylim([-0.4;0.6]);xlim([20;100]);
  
        dataBreak_nodes(ITERbn,1)=t;
        dataBreak_nodes(ITERbn,2:length(Break_nodes)+1)=Break_nodes;
        dataCrestBreak(ITERbn,1)=t;
        dataCrestBreak(ITERbn,2:length(CrestBreakNew.indexNow)+1)=CrestBreakNew.indexNow;
        ITERbn=ITERbn+1;
    end
    
     CrestBreakNew=[];
end

if t>t_init+ITERbdt*param.dt   
        ITERbdt=ITERbdt+1;
end



