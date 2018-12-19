function [heat]=heatindex(tmpf,relh);

% computes the NWS heat indices using formulas from Steadman (1979)
% heat is heat index in degrees F
% heat is prescribed values of
% 0 is there is missing data
% 1 when temp<40F
% 2 when temp>87F and rh>85F (no-analog)
% 3 when temperature and RH exceed values along the diagonal (no-analog)
% 4 when temp>112F and RH>13% (no-analog)
% all other values are recast to unsigned integers for computation
% efficiency


% /************************************************************************
%  * pd_heat								*
%  *									*
%  * This subroutine computes HEAT, the Southern Region/CPC Rothfusz heat	*
%  * index, from TMPF and RELH.  The output will be calculated in degrees *
%  * Fahrenheit.								*
%  *									*
%  *	Source:  NWS Southern Region SSD Technical Attachment SR 90-23	*
%  * 		 7/1/90.  Heat Index was originally known as the	*
%  *		 apparent temperature index (Steadman, JAM, July, 1979).*
%  *									*
%  * The Rothfusz regression is optimal for TMPF > ~80 and RELH > ~40%.	*
%  * This code applies a simple heat index formula and then resorts to	*
%  * the Rothfusz regression only if the simple heat index exceeds 80,	*
%  * implying temperatures near, but slightly below 80.  To make the	*
%  * simple calculation continuous with the values obtained from the	*
%  * Rothfusz regression, the simple result is averaged with TMPF in	*
%  * computing the simple heat index value.				*
%  *									*
%  * This code includes adjustments made by the CPC for low RELH at high	*
%  * TMPF and high RELH for TMPF in the mid 80's.				*
%  *									*
%  * pd_heat ( tmpf, relh, np, heat, iret )				*
%  *									*
%  * Input parameters:							*
%  *	*tmpf		const float	Temperature in Fahrenheit	*
%  *	*relh		const float	Relative humidity in percent	*
%  *	*np		const int	Number of points		*
%  *									*
%  * Output parameters:							*
%  *	*heat		float		Heat Index in Fahrenheit	*
%  *	*iret		int		Return code			*
%  *					  0 = normal return		*
%  **									*
%  * Log:									*
%  * S. Jacobs/NCEP	 3/02	Initial coding				*
%  * K. Brill/HPC		11/02	Remove SUBFLG input logical array	*
%  * K. Brill/HPC		 1/03	Fix discontinuity around 77 F		*
%  * R. Tian/SAIC		 9/05	Translated from FORTRAN			*
%  J. Abatzoglou 12/18 Translated to MATLAB
%  ************************************************************************/

%If the temperature is less than 40 degrees, then set the heat index to the temperature.
relh=relh';tmpf=tmpf';
heat=(zeros(size(tmpf)));
f=find(tmpf<=40);heat(f)=tmpf(f);

% if temperature is above 40 then calculate heat indices
f=find(tmpf>40);

%Compute a simple heat index. If the value is less than 80 degrees, use it.

heat(f) = 61.0 + (tmpf(f)-68.0) * 1.2 + relh(f) * .094;
heat(f) = .5 .* ( tmpf(f) + heat(f) );

%For values above 80F compute full regression value of the heat index.

f=find(heat>=80);

t2 = tmpf(f).*tmpf(f);
r2 = relh(f).*relh(f);
heat(f) = -42.379+  2.04901523 * tmpf(f)+ 10.14333127 * relh(f) -  0.22475541 * tmpf(f) .* relh(f)-  0.00683783 * t2 -  0.05481717 * r2+  0.00122874 * t2 .* relh(f)+  0.00085282 * tmpf(f) .*r2 -  0.00000199 * t2 .* r2;

% Adjust for high regression at low RH for temps above 80 degrees F.
f2=find(relh <=13 & tmpf>=80 & tmpf<=112);
adj1 = ( 13.0 - relh(f2) ) / 4.0;
tval = 17.0 - abs ( tmpf(f2) - 95.0 );
adj2 = sqrt ( tval / 17.0 );
adj  = adj1 .* adj2;
heat(f2) = heat(f2)-adj;

%  Adjust for low regression at high RH and temps in the mid 80s.
f2=find(relh > 85.0 & tmpf >= 80.0 & tmpf <= 87.0 );
adj1 = ( ( relh(f2) - 85.0 ) / 10.0 );
adj2 = ( ( 87.0 - tmpf(f2) ) / 5.0 );
adj  = adj1 .* adj2;
heat(f2) = heat(f2)+adj;

% in=3 for invalid data along the diagonal
f=find(tmpf>102 & relh>40);heat(f)=3;
f=find(tmpf>101 & relh>50);heat(f)=3;
f=find(tmpf>100 & relh>55);heat(f)=3;
f=find(tmpf>98 & relh>60);heat(f)=3;
f=find(tmpf>97 & relh>65);heat(f)=3;
f=find(tmpf>96 & relh>70);heat(f)=3;
f=find(tmpf>95 & relh>75);heat(f)=3;
f=find(tmpf>92 & relh>80);heat(f)=3;

% in=2 for muggy+hot conditions temp>87 and rh>85
f=find(tmpf>87 & relh>=85);heat(f)=2;

% in=4 for really hot temperatures with rh>13
f=find(tmpf>112 & relh>13);heat(f)=2;

% in=1 for cold temperatures temp<40
f=find(tmpf<40);heat(f)=1;

heat=uint8(heat);