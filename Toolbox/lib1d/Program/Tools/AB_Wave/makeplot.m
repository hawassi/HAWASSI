function makeplot(input1,input2,x,time,eta,savename,dynmodel,folder_root)
% input1 : profile or signal
% input2 : tsnap or Xbouy
if  strcmp(input1,'xprofile')
    tsnap = input2;
    if tsnap > max(time)  ,  tsnap = max(time);
    end
                    MTC     = max(eta);  MTT     = min(eta);
            figure(1000)
                plot(x, eta(closest(time,tsnap),:),'b',...
                        x, MTC, 'k',x, MTT, 'c',x,0,'k--');
                xlabel('x [m]') ; ylabel('\eta [m]');
%                 grid on;
            title([savename,', ',dynmodel,', MTA & profile @ time =',num2str(tsnap),'[s]']);                   
            saveas(gcf,[folder_root,'\Output\Profile_t',num2str(tsnap),'.fig'])
                 
elseif strcmp(input1,'tsignal')
    Xbuoy =input2;
    if Xbuoy > max(x)  ,  Xbuoy = max(x);
    end            
                    MSC     = max(eta,[],2)';   MST  = min(eta,[],2)';
                    buoyind     = closest(x,Xbuoy);
                    etabuoy = eta(:,buoyind);
            figure(1001)
                plot(time, etabuoy,'r',...
                        time,MSC,'k',time,MST,'c',time,0,'k--');
                xlabel('time [s]') ;   ylabel('\eta [m]');
%                 grid on;
            title([savename,', ',dynmodel,', MTA & time signal @ position =',num2str(Xbuoy),'[m]']);
%                 axis([300 13000 -3 5]);                   
            saveas(gcf,[folder_root,'\Output\Signal_x',num2str(Xbuoy),'.fig'])
                
end
end