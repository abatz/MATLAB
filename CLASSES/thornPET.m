function [PET]=thornPET(T,lat);
%function [PET]=thornPET(T,lat);
% function estimates monthly potential evapotranspiration using the Thorthwaitte method
% requires monthly data (12-months) and latitude of observation (to deduce day length)
% function returns monthly PET in units of mm

f=find(T<0);T(f)=0;
I=zeros(12,length(lat));
for i=1:12
 I(i,:)=(T(i,:)/5).^1.514;
end
I=sum(I,1);
alpha=6.75e-7*I.^3-7.71e-5*I.^2+1.792e-2*I+.49239;


d2=[31 59 90 120 151 181 212 243 273 304 334 365];
d1=[1 32 60 91 121 152 182 213 244 274 305 335];
daysinmonth=d2-d1+1;
for month=1:12
for i=1:daysinmonth(month)
 L(i,:)=calcDay(i-1+d1(month),lat);
end
L=nmean(L,1);
L=L/12;
N=daysinmonth(month)/30;
TERM3=(10*T(month,:)./I).^alpha;
PET(month,:)=1.6.*L.*N.*TERM3;
end
PET=PET*10;


function [daylit]=calcDay(j_date,latdec);

   phi = 0;
   xfact = 0;
   decl = 0;
   tla = 0;
 
f=find(j_date>365);
j_date(f)=365;
phi=latdec*.01745;
decl=.41008*sin((j_date-82)*.01745);
daylit=24*(1-acos(tan(phi)*tan(decl))/3.14159);daylit=real(daylit);

