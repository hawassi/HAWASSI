function IdCheckAddCrest=funBr_check_addcrest(kk,CrestBreak,CrestBreakPrev)
NCBprev=length(CrestBreakPrev.index);
IndexNewCrest=CrestBreak.index(kk);
IdCheckAddCrest=1;
for ll=1:NCBprev
    if  any(CrestBreakPrev.ID(ll).nodes==IndexNewCrest)
        IdCheckAddCrest=0; %restrict to different wave
        break;
    end
end