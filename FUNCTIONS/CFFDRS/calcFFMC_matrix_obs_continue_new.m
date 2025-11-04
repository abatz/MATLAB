function [ffmc]=calcFFMC(pt,tmx,rh,wind,season);
% #############################################################################
%   # Description: Fine Fuel Moisture Code Calculation. All code
%   #              is based on a C code library that was written by Canadian
%   #              Forest Service Employees, which was originally based on
%   #              the Fortran code listed in the reference below. All equations
%   #              in this code refer to that document.
%   #
%   #              Equations and FORTRAN program for the Canadian Forest Fire 
%   #              Weather Index System. 1985. Van Wagner, C.E.; Pickett, T.L. 
%   #              Canadian Forestry Service, Petawawa National Forestry 
%   #              Institute, Chalk River, Ontario. Forestry Technical Report 33. 
%   #              18 p.
%   #
%   #              Additional reference on FWI system
%   #
%   #              Development and structure of the Canadian Forest Fire Weather 
%   #              Index System. 1987. Van Wagner, C.E. Canadian Forestry Service,
%   #              Headquarters, Ottawa. Forestry Technical Report 35. 35 p.
%   #  
%   #
%   # Args: ffmc_yda:   The Fine Fuel Moisture Code from previous iteration
%   #           temp:   Temperature (centigrade)
%   #             rh:   Relative Humidity (%)
%   #           prec:   Precipitation (mm)
%   #             ws:   Wind speed (km/h)
%   #       
%   #
%   # Returns: A single ffmc value
%   #
%   #############################################################################

wind=wind*3.6;
s=size(pt);
mr=zeros(size(pt));
pt=shiftdim(reshape(shiftdim(reshape(pt,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
tmx=shiftdim(reshape(shiftdim(reshape(tmx,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
rh=shiftdim(reshape(shiftdim(reshape(rh,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
wind=shiftdim(reshape(shiftdim(reshape(wind,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
season=shiftdim(reshape(shiftdim(reshape(season,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);

FFMCY=85*ones(size(pt,1),1);lastseasony=FFMCY*0;
for yr=1:size(pt,3);  % loop over yr
for i=1:size(pt,2);  % loop over days
%   #Eq. 1
  wmo = 147.2 * (101 - FFMCY)./(59.5 + FFMCY);
%   #Eq. 2 Rain reduction to allow for loss in 
%   #  overhead canopy
  ra=max(0,pt(:,i,yr)-0.5);
%ra=pt(:,i,yr)-0.5;
%   #Eqs. 3a & 3b
  wmonew=wmo;
  f=find(pt(:,i,yr)>0.5 & wmo>150);
  wmonew(f)=wmo(f)+0.0015 * (wmo(f) - 150) .* (wmo(f) - 150) .* sqrt(ra(f)) + 42.5 .* ra(f) .* exp(-100 ./ (251 - wmo(f))) .* (1 - exp(-6.93 ./ ra(f)));
  f=find(pt(:,i,yr)>0.5 & wmo<=150);
  wmonew(f)=wmo(f) + 42.5 * ra(f) .* exp(-100 ./ (251 - wmo(f))) .*(1 - exp(-6.93 ./ ra(f)));
  wmonew=reshape(wmonew,size(wmo));
%   #The real moisture content of pine litter ranges up to about 250 percent,
%   # so we cap it at 250  
  wmo=min(wmonew,250);
%   #Eq. 4 Equilibrium moisture content from drying
  ed = 0.942 * (rh(:,i,yr).^0.679) + (11 * exp((rh(:,i,yr) - 100) / 10)) + 0.18 * (21.1 - tmx(:,i,yr)) .* (1 - 1 ./ exp(rh(:,i,yr) * 0.115));
%   #Eq. 5 Equilibrium moisture content from wetting
  ew = 0.618 * (rh(:,i,yr).^0.753) + (10 * exp((rh(:,i,yr) - 100) / 10)) + 0.18 * (21.1 - tmx(:,i,yr)) .* (1 - 1 ./ exp(rh(:,i,yr) * 0.115));
%   #Eq. 6a (ko) Log drying rate at the normal
%   #  termperature of 21.1 C
  f=find(wmo<ed & wmo<ew);
  z=zeros(size(wmo));

  z(f)= 0.424 .* (1 - (((100 - rh(f,i,yr)) ./ 100).^1.7)) + 0.0694 .* sqrt(wind(f,i,yr)) .* (1 - ((100 - rh(f,i,yr)) ./ 100).^8);
%   #Eq. 6b Affect of temperature on  drying rate
  x =z * 0.581 .* exp(0.0365 .* tmx(:,i,yr));
%   #Eq. 8
  wm=wmo;
  f=find(wmo<ed & wmo<ew);
  wm(f)=ew(f) - (ew(f) - wmo(f))./(10.^x(f));
  
%   #Eq. 7a (ko) Log wetting rate at the normal
%   #  termperature of 21.1 C    
  f=find(wmo>ed);
  z(f)= 0.424 .* (1 - (rh(f,i,yr)/100).^1.7) + 0.0694 .* sqrt(wind(f,i,yr)) .* (1 - (rh(f,i,yr)/100).^8);
%   #Eq. 7b Affect of temperature on  wetting rate
  x = z .* 0.581 .* exp(0.0365 .* tmx(:,i,yr));
%   #Eq. 9
  f=find(wmo>ed);
  wm(f)=ed(f)+(wmo(f)-ed(f))./(10.^x(f));

%   #Eq. 10 Final ffmc calculation
  ffmc1 = (59.5 * (250 - wm))./(147.2 + wm);
%   #Constraints
  ffmc1=min(ffmc1,101);
  ffmc1=max(ffmc1,0);
  ffmc(:,i,yr)=ffmc1;FFMCY=ffmc1;
  f=find(season(:,i,yr)==0);
  ffmc(f,i,yr)=50;FFMCY(f)=50;
  f=find(lastseasony==0 & season(:,i,yr)==1);FFMCY(f)=85; % standard season start DC no-overwintering
  lastseasony=season(:,i,yr);
end
end
ffmc=real(ffmc);
