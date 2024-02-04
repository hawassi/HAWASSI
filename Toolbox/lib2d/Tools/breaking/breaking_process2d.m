%Version 1 Sept 2016
function [B,AllBreakNodes] = breaking_process2d(eta,gradphi,dom,param,tnow,g,t_init)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Passing parameters
U       = sqrt(gradphi.x.^2+gradphi.y.^2);
B       = zeros(size(eta));

global flagbr iterNprev  CrestBreakPrev ...
    ITERbn dataBreak_nodes ITERbdt   dataCrestBreak
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Kinematic criteria U>C and find all peaks
[CrestBreak] = Kinematic_breaking_criterion2d(tnow,eta,gradphi,U,dom,param,g);

if ~isempty(CrestBreak)&&flagbr==0
    flagbr=1;
    CrestBreakNew=CrestBreak;
elseif isempty(CrestBreak) &&flagbr==2
    CrestBreakNew=CrestBreakPrev;
elseif ~isempty(CrestBreak)&& flagbr==2
    %%
    NCB=length(CrestBreak.index) ;
    NCBprev=length(CrestBreakPrev.index);
    CrestBreakNew=CrestBreakPrev;
    
    add=1;
    for kk=1:NCB 
        IdCheckAddCrest=funBr_check_addcrest(kk,CrestBreak,CrestBreakPrev);

        if IdCheckAddCrest ==1 % It will be saved if Crest break candidates are not in the previous breaking nodes 
            CrestBreakNew.index(NCBprev+add)   =CrestBreak.index(kk);%|
            CrestBreakNew.indexNow(NCBprev+add)=CrestBreak.index(kk);%|same for the first time
            CrestBreakNew.ID(NCBprev+add).nodes=CrestBreak.index(kk);
            CrestBreakNew.Ucrest(NCBprev+add)  =CrestBreak.Ucrest(kk);
            CrestBreakNew.position.x(NCBprev+add)=CrestBreak.position.x(kk);
            CrestBreakNew.position.y(NCBprev+add)=CrestBreak.position.y(kk);
            CrestBreakNew.time(NCBprev+add)    =CrestBreak.time(kk);
            CrestBreakNew.DirProp.theta(NCBprev+add) =CrestBreak.DirProp.theta(kk);
            CrestBreakNew.DirProp.kwadran(NCBprev+add) =CrestBreak.DirProp.kwadran(kk);
            add=add+1;
        end
    end
end

if (flagbr==1 || flagbr==2)
   NCBnow=length(CrestBreakNew.index);  
    CrestBreakPrev=[];
    iterNprev=1;
    IDTerminate=zeros(NCBnow,1);
     
    AllBreakNodes=[]; iterbn=1;
    for lx=1:NCBnow  %% looping in index of peaks
        if flagbr==1
            CrestBreakNew.indexNow(lx)=CrestBreakNew.index(lx);
        else
            indxMaxNow=funBr_find_peak_loc(param,eta,dom,CrestBreakNew,lx);
            if isempty(indxMaxNow)
               IDTerminate(lx)=1;
               continue; 
            end
            CrestBreakNew.indexNow(lx)=indxMaxNow;
            CrestBreakNew.position.x(lx)=dom.XX(indxMaxNow);
            CrestBreakNew.position.y(lx)=dom.YY(indxMaxNow);
            DirPropTheta=atan(gradphi.y(indxMaxNow)/gradphi.x(indxMaxNow));
            kwadran=funBr_check_kwadran_dir(gradphi.x(indxMaxNow),gradphi.y(indxMaxNow));
            CrestBreakNew.DirProp.theta(lx)=abs(DirPropTheta);
            CrestBreakNew.DirProp.kwadran(lx)=kwadran;
        end
   
        tbreak=CrestBreakNew.time(lx);
        U_I=CrestBreakNew.Ucrest(lx);
        U_F=param.TC.*U_I;
        if  tnow-tbreak>=param.Tstar
            U_star=U_F;
        elseif (tnow-tbreak >=0) && (tnow-tbreak <param.Tstar)
            U_star=U_I+((tnow-tbreak)/param.Tstar).*(U_F-U_I);
        else
            U_star=U_I;
        end
        
       [B, Break_nodes]=funBr_BreakNodes_char(param,U,U_star,U_F,eta,dom,CrestBreakNew,B,lx);
    
        if all(U(Break_nodes)<U_star) %&& t-tbreak>0
            IDTerminate(lx)=1;
        else
            IDTerminate(lx)=0;
            CrestBreakPrev.ID(iterNprev).nodes=Break_nodes;
            CrestBreakPrev.index(iterNprev)=CrestBreakNew.index(lx);
            CrestBreakPrev.indexNow(iterNprev)=CrestBreakNew.indexNow(lx);
            CrestBreakPrev.Ucrest(iterNprev)=CrestBreakNew.Ucrest(lx);
            CrestBreakPrev.position.x(iterNprev)=CrestBreakNew.position.x(lx);
            CrestBreakPrev.position.y(iterNprev)=CrestBreakNew.position.y(lx);
            CrestBreakPrev.time(iterNprev)=CrestBreakNew.time(lx);
            CrestBreakPrev.DirProp.theta(iterNprev)=CrestBreakNew.DirProp.theta(lx);
            CrestBreakPrev.DirProp.kwadran(iterNprev)=CrestBreakNew.DirProp.kwadran(lx);
            iterNprev=iterNprev+1;
        end
        Nnodes=length(Break_nodes);
        AllBreakNodes(iterbn:iterbn+Nnodes-1)=Break_nodes;
        iterbn=iterbn+Nnodes;
    end
%     global iterdtcheck
%     if tnow>86 && tnow<90     
% %             subplot(3,1,3)
%         if mod( iterdtcheck,2)
%             surf(dom.XX,dom.YY,eta,'Edgecolor','none')
%             xlabel('x');ylabel('y');view([30,60]);
%             hold on;
%             quiver(dom.XX(AllBreakNodes),dom.YY(AllBreakNodes),gradphi.x(AllBreakNodes),gradphi.y(AllBreakNodes),5);
%             plot_properties;
%             hold on;
%             plot3(dom.XX(AllBreakNodes),dom.YY(AllBreakNodes),eta(AllBreakNodes),'ow')
%             xlabel('x');ylabel('y');          
%             hold off
%             pause(0.001)
%         end
%     end
    
    if all(IDTerminate==1)
        flagbr=0;
    else
        flagbr=2;
    end
    
    if tnow>t_init+ITERbdt*param.dt  
        dataBreak_nodes(ITERbn,1)=tnow;
        dataBreak_nodes(ITERbn,2:length(AllBreakNodes)+1)=AllBreakNodes;
        dataCrestBreak(ITERbn).t=tnow;
        dataCrestBreak(ITERbn).X=CrestBreakNew.position.x;
        dataCrestBreak(ITERbn).Y=CrestBreakNew.position.y;
        ITERbn=ITERbn+1;
    end
end

% 
 if tnow>t_init+ITERbdt*param.dt   
     ITERbdt=ITERbdt+1;
 end

