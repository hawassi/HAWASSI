% invOm, version March 2012
% % % % Software developed by 
% % % % LabMath-Indonesia & University of Twente
function Ddeep = invUpCgDeep(v,D0,deltaD,UgCgrat)
Ddeep=D0;
k=invOmExact(v,Ddeep);
while UgExact(k,Ddeep)/UpExact(k,Ddeep)>UgCgrat
    k=invOmExact(v,Ddeep);
    Ddeep=Ddeep+deltaD;
end