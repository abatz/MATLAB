function [ISI]=isi_calc(ffmc,vs);
vs=vs*3.6;  %m/s to kph
%   ==========================================================================  isi_calc
%  |  Subroutine:  isi_calc                                                   | isi_calc
%  |  Author:      ???/CFS:YYYYMMDD                                           | isi_calc
%  |  Modified:    Modifications listed start with the most recent:           | isi_calc
%  |               RJC/CFS:20060404  Added Imperial option, label and comments| isi_calc
%  |  Purpose:     Calculate Initial Spread Index (ISI) for the CFFWIS.       | isi_calc
%  |  Args:        $Ffmc -- CFFWIS Fine Fuel Moisture Code                    | isi_calc
%  |               $Wind -- Wind speed (km/h)                                 | isi_calc
%  |               $Imperial -- Wind speed in Imperial (1) or SI (0) units    | isi_calc
%  |  Globals:     none                                                       | isi_calc
%  |  Returns:     $R -- isi                                                  | isi_calc
%   --------------------------------------------------------------------------  isi_calc
%
 
    fW = exp(0.05039 .* vs);                                          
    m = 147.2 * ((101.0 - ffmc)./(59.5 + ffmc));                       
    FF = 91.9 * exp(-.1386 .* m) .* ( 1 + (reshape(m(:).^5.31,size(m)))/4.93E7);           
    R = 0.208 .* fW .* FF;                                           
    f=find(R<0);R(f)=0;
       ISI=reshape(R,size(ffmc));
