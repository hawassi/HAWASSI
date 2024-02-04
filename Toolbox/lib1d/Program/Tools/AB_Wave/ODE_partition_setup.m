auto_partition=options.partition.auto;
tintval=influx.gen.timesig;
if strcmp(auto_partition,'on')
    Nx=length(x);Nt=length(tintval);
    Matrix_Size=Nx*Nt*2; %eta and u;
    factdata=1000000;
    [c,max_matrix_CPU]=computer;
    if Matrix_Size> max_matrix_CPU/factdata
        Npartition=floor(Matrix_Size./(max_matrix_CPU/factdata));
    else
        Npartition=1;
    end
else
    Npartition=options.partition.N;
end
  par.ode_Npartition=Npartition;