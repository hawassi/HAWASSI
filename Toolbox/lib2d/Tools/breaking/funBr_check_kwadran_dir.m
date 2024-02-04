function kwadran=funBr_check_kwadran_dir(ux,uy)
if ux>=0 &&uy>=0
    kwadran=1;
elseif ux<0 &&uy>=0
    kwadran=2;
elseif ux<0 && uy<0
    kwadran=3;
else
    kwadran=4;
end