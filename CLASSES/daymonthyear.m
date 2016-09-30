function [d,m,y]=daymonthyear(doy,year);
d2=[31 59 90 120 151 181 212 243 273 304 334 365];
d1=[1 32 60 91 121 152 182 213 244 274 305 335];

if doy>60 & mod(year,4)==0
   d1(3:12)=d1(3:12)+1;
   d2(2:12)=d2(2:12)+1;
end
      
for i=1:length(doy)
    %if mod(year(i),4)==0
        f=find(d1<=doy(i) & d2>=doy(i));
if ~isempty(f)
        m(i)=f;d(i)=1+doy(i)-d1(f);y(i)=year(i);
else
 d=NaN;m=NaN;y=NaN;
end
end

 
