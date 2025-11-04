function [BUI]=bui_calc(dmc,dc);

%
%   ==========================================================================  bui_calc
%  |  Subroutine:  bui_calc                                                   | bui_calc
%  |  Author:      ???/CFS:YYYYMMDD                                           | bui_calc
%  |  Modified:    Modifications listed start with the most recent:           | bui_calc
%  |               RJC/CFS:20060404  Added label and comments                 | bui_calc
%  |  Purpose:     Calculate BuildUp Index (bui) for the CFFWIS.              | bui_calc
%  |  Args:        $dmc -- CFFWIS Duff Moisture Code                          | bui_calc
%  |               $dc -- CFFWIS Drought Code                                 | bui_calc
%  |  Globals:     none                                                       | bui_calc
%  |  Returns:     $bui -- bui                                                | bui_calc
%   --------------------------------------------------------------------------  bui_calc
%
%
% ............................................................................|
%    BUI calculation ... numbers appearing as comments are equation numbers in
%    Van Wagner (1985) Equations and FORTRAN Program for the Canadian Forest
%    Fire Weather Index System:
%	
dd=dmc+0.4*dc;
BUI=NaN*ones(size(dd));
f1=find(dd~=0);f2=find(dd==0);
BUI(f2)=0;
f3=find(dmc(f1)<=0.4*dc(f1));
BUI(f1(f3)) = (0.8 * dmc(f1(f3)) .* dc(f1(f3)))./(dd(f1(f3)));             % 27a
f4=find(dmc(f1)>0.4*dc(f1));
s=size(dmc);
dmc=dmc(f1(f4));dc=dc(f1(f4));dd=dd(f1(f4));
BUI(f1(f4)) = dmc - (1 - (0.8 .* dc ./dd)).*(0.92 + reshape((.0114.*dmc).^1.7,[length(f4) 1]));
BUI=reshape(BUI,s);

