function [AET,DEF,RUNOFF]=simplehydromodel(PPT,PET,AWC);
%function [AET,DEF,RUNOFF]=simplehydromodel(PPT,PET,AWC);
% this function solves the monthly water balance code using a very simple
% one bucket model and assumes that all water is liquid.
% PPT and PET should be in mm 
% optional argument: AWC = soil water holding capacity and should be in mm, default is 150mm
% returns monthly AET (actual evapotranspiration), DEF (unmet demand, or deficit) and RUNOFF
% all values returned are in units of mm
% as a default the initial month starts with 0 soil moisture and is spun up over a year before
% retaining actual values ONLY for the second year


AET=NaN*ones(size(PPT));DEF=AET;RUNOFF=AET;
if nargin==2
AWC=150;
end
soil=0;
for j=1:2
for i=1:length(PPT);
deltasoil=PPT(i)-PET(i);
if soil+deltasoil<0
 AET(i)=soil+PPT(i);
 DEF(i)=-(soil+deltasoil);
 RUNOFF(i)=0;
 soil=0;
else
 AET(i)=PET(i);
 DEF(i)=0;
 if soil+deltasoil>AWC
  RUNOFF(i)=soil+deltasoil-AWC;
  soil=AWC;
 else
  soil=soil+deltasoil;
   RUNOFF(i)=0;
 end
end
end;end

