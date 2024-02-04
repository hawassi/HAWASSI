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

function main(GUISAmain)
proj.dir=GUISAmain.proj.projdir;
proj.name=GUISAmain.proj.name;
proj.workdir=GUISAmain.proj.workdir;

ship.shape=GUISAmain.ship.shape;
ship.shapedat=GUISAmain.ship.shapedat;
ship.width=GUISAmain.ship.width;
ship.draft=GUISAmain.ship.draft;
ship.cog=GUISAmain.ship.cog;
ship.Nx=GUISAmain.ship.Nx;
ship.visc_sway=GUISAmain.ship.visc_sway;
ship.visc_roll=GUISAmain.ship.visc_roll;
ship.GM=GUISAmain.ship.GM;
ship.gyradius=GUISAmain.ship.gyradius;
ship.mooring.check=GUISAmain.ship.mooring.check;
ship.mooring.Tn=GUISAmain.ship.mooring.Tn;

wave.inputVar=GUISAmain.wave.freqVar;
wave.input=GUISAmain.wave.freq;
wave.amplVar=GUISAmain.wave.amplVar;
wave.ampl=GUISAmain.wave.ampl;

bottom.depth=GUISAmain.bottom.depth;
bottom.slope=GUISAmain.bottom.slope;

calc.Id=GUISAmain.calc.Id;% 0 diffraction, 1 raditaion, 2 RAO
calc.evmodes=GUISAmain.calc.evmodes;
calc.valdata=GUISAmain.data;
calc.normVartype=GUISAmain.normVartype;

calc.difraddata=GUISAmain.dataDiffRad;

%%%
if ship.shape==2
 ship.width=2*ship.draft;
end
[par]=funSA_preparation(wave,ship,bottom,calc);
assignin('base','par',par);

if ship.shape==4
  ship.width=par.S4B;
  ship.draft=par.S4draft;
end
    
if calc.Id==1||calc.Id==2
    [OutdiffRadd]=funSA_diffraction_or_radiation(par,ship,wave,bottom,calc);
elseif calc.Id==3
    [OutRAO]=funSA_RAO(proj,par,ship,wave,bottom,calc);
end 

if calc.Id==1
    funSA_ploting_diffraction;
elseif  calc.Id==2
    funSA_ploting_radiation;
elseif calc.Id==3
    funSA_ploting_RAO;
end
end
