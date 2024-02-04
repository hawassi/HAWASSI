%% Properties of dispersion: Exact and Used  
function DispModel(k,depth,Omd,Ugd,fundispersion,sf_savename)
%     depth = par.depth;
%     Om      = str2func(fun.dispersion);      
%     Omd     = @(k)Om(k,depth);
%     Ug      = str2func(fun.groupvel);
%     Ugd     = @(k)Ug(k,depth);
    kk  = k(1:floor(length(k)/2)-1);
    dk = k(3)-k(2);
figure(102)
    subplot(2,1,1)
        plot(kk,OmExact(kk,depth),'b',kk,Omd(kk),'r--')
        axis([0 max(k)-10*dk 0 Omd(max(k))]);
        grid on;
        title(['Exact Dispersion (blue) and ',fundispersion,' (red-);'...
                    ' depth = ',num2str(depth),' [m]'])
        ylabel('frequency [radian/s]');
    subplot(2,1,2)
        plot(kk,UgExact(kk,depth),'b',kk,Ugd(kk),'r--');
        axis([0 (max(k)-10*dk) 0 Ugd(k(1))]);
        grid on;    
        title('Group velocity')
        xlabel('wave number [rad/m]') ;   ylabel('group vel [m/s]');           
        saveas(gcf,[sf_savename,'DispersionComp','.fig'])
end        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


