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

if strcmpi(input.wave.option,'yes')
Source_hat=0; 
for I=1:influx.input.N
    INsig=influx.gen(I).Signal;
    Gam=influx.gen(I).Gam;
    Indirection=model.influx.direction(I);
        INsig_t0      =ppval(INsig,time);
    if strcmp(influx.spatial.line(I).Orientation,'Vertical')
        INsig_t      =repmat(INsig_t0,1,dom.Nx);
    else
        INsig_t      =repmat(INsig_t0,1,dom.Ny)';
    end
    
    if ~strcmp(Indirection, 'Bi')
    SkewGam      = influx.gen(I).SkewGam;
    INsig_skew   = influx.gen(I).SkewSignal;
    SkewINsig_t  = ppval(INsig_skew,time);
    if strcmp(influx.spatial.line(I).Orientation,'Vertical')
        SkewINsig_t      =repmat(SkewINsig_t,1,dom.Nx);
    else
        SkewINsig_t      =repmat(SkewINsig_t,1,dom.Ny)';
    end
   
    alpha        = 1;  
    if strcmp(Indirection, 'UniR')||strcmp(Indirection, 'UniT')
    beta=-1;
    elseif strcmp(Indirection, 'UniL')||strcmp(Indirection, 'UniB')
    beta=1;    
    end     
     
    SkewSource_hat= fft2(beta.*SkewINsig_t.*SkewGam);
    else
    SkewSource_hat  = 0; alpha        = 2;
    end
     Source         =alpha.*INsig_t.*Gam; 
     Source_hat_I   = fft2(Source)+SkewSource_hat;
     Source_hat     =Source_hat+Source_hat_I;
%      SourceSkew=funC_ifft2(SkewSource_hat);
end
else
    Source_hat=0;
end