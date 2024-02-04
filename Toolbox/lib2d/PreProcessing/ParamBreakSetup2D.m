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

%% parameters of breaking
%these 4 parameters below have been declared in MainUser.m
%par.break.KBC       = GUImain.Break_param(1); % Kin.Br.Cond.: U/C in [0.7:1]
%par.break.TC        = GUImain.Break_param(2); % Terminal Condition uF=TC*uI; TC in [0.2:0.3]
%par.break.Tchar     = GUImain.Break_param(3); % Characteristic Time T*=Tchar*Tp; Tchar in [0.1:0.5]; Tp->peak period
%par.break.delb      = 1.2              ; % Mixing length Coef
%Another parameters are below
if strcmp(model.breaking.check,'Yes')
    parBreak.KBC =model.breaking.KBC;
    parBreak.TC  =model.breaking.TC;
    parBreak.delb=model.breaking.delb;
    
    if strcmpi(input.wave.option,'Yes')
        parBreak.MaxEtaInit=max(influx.par.Hs)/2*0.9;               % as reference to find peaks (Max)
        parBreak.Tstar = model.breaking.Tchar*min(influx.par.T_p);  %characteristic time breaking [Tp/10;Tp/2]
        parBreak.Rsearch   =2*pi./max(influx.par.k_p)./4;
    elseif ivp.type>1
        parBreak.MaxEtaInit=max(max(ivp.eta))/2;               % as reference to find peaks (Max)
        parBreak.Tstar = 5*sqrt(ivp.par.meandepth*g);             %based on Kennedy
        parBreak.Rsearch   =2*pi./ivp.par.k_p./4;
    elseif input.bdyassim.option==1
        parBreak.MaxEtaInit=max(bdyassim.par.Hs)/2*0.9;              % as reference to find peaks (Max)
        parBreak.Tstar = 5*sqrt(bdyassim.par.meandepth*g);             %based on Kennedy
        parBreak.Rsearch   =2*pi./bdyassim.par.k_p./4;
    end
    parBreak.dt        =timeSimul.interval(2)-timeSimul.interval(1);  %for saving data each dt
    parBreak.char      =dom.cfSA;
    
    parBreak.char(dom.cfSA<1)=0;
    parBreak.char(influx.ChiAdj<1)=0;
    
    if  strcmpi(input.wall.option,'Yes')
        parBreak.char(dom.wall.char<1)=0;
    end
    
    if  input.body.option==1
        parBreak.char(body.par.chiB(end).char>0)=0;
    end
    
else
    parBreak=[];
end

% size(parBreak.char)
% figure(11);
% set(gcf,'Renderer','zbuffer'); %due to graphics driver
% surf(dom.XX,dom.YY,parBreak.char,'edgecolor','none')
%  view(2)