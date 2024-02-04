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

auto_partition=options.partition.auto;
tintval=timeSimul.interval;
if strcmp(auto_partition,'on')
    Nx=dom.Nx;Ny=dom.Ny;Nt=length(tintval);
    Matrix_Size=Nx*Ny*Nt*3; %eta and u;
    factdata=2000000;
    [c,max_matrix_CPU]=computer;
    if Matrix_Size> max_matrix_CPU/factdata;
        Npartition=floor(Matrix_Size./(max_matrix_CPU/factdata));
    else
        Npartition=1;
    end
else
    Npartition=options.partition.N;
end
  par.ode.Npartition=Npartition;