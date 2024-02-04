function[Ndom]=funC_Ndom(NdxS)
    Ndom=(NdxS-1)/2+2; % 10; % number of domain 16
    if mod(Ndom,2)~=0
        Ndom=Ndom+1;
    end
end