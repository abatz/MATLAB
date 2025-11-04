    function [FWI]=calcFWI(ISI,BUI);
%
%   ==========================================================================  fwi_calc
%  |  Subroutine:  fwi_calc                                                   | fwi_calc
%  |  Author:      ???/CFS:YYYYMMDD                                           | fwi_calc
%  |  Modified:    Modifications listed start with the most recent:           | fwi_calc
%  |               RJC/CFS:20060404  Added label and comments                 | fwi_calc
%  |  Purpose:     Calculate Fire Weather Index (FWI) for the CFFWIS.         | fwi_calc
%  |  Args:        $isi -- CFFWIS Initial Spread Index                        | fwi_calc
%  |               $bui -- CFFWIS BuildUp Index                               | fwi_calc
%  |  Globals:     none                                                       | fwi_calc
%  |  Returns:     $S -- fwi                                                  | fwi_calc
%   --------------------------------------------------------------------------  fwi_calc
%
% ............................................................................|
%    FWI calculation ... numbers appearing as comments are equation numbers in
%    Van Wagner (1985) Equations and FORTRAN Program for the Canadian Forest
%    Fire Weather Index System:
%
f1=find(BUI<=80);
fD(f1)=0.626.*BUI(f1).^.809+2;
f2=find(BUI>80);
fD(f2)= 1000.0./(25 + (108.64 * exp(-0.023.* BUI(f2))));                % 28b
f3=find(isnan(BUI));fD(f3)=NaN;
B = 0.1 .* ISI(:) .*fD(:);                                                % 29
f1=find(B>1);f2=find(B<=1);f3=find(isnan(B));
FWI(f1)= exp(2.72 .* ((0.434*log(B(f1))).^0.647));                        % 30a
FWI(f2)=B(f2);
FWI(f3)=NaN;
FWI=reshape(FWI,size(ISI));
