function [dmc]=calcDMC(pt,tmx,rh,lat,season);


s=size(pt);
lastseasony=zeros(size(pt,1),1);
if s(1)==1|s(2)==1

k1=repmat((1:365)',[1 size(tmx,1)]);
k2=repmat(lat,[365 1]);
day=real(day_length(k1,k2));
else
day=day_length(repmat((1:365)',[1 size(tmx,1) size(tmx,2)]),shiftdim(repmat(lat,[1 size(tmx,2) 365]),2));
end
pt=shiftdim(reshape(shiftdim(reshape(pt,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
tmx=shiftdim(reshape(shiftdim(reshape(tmx,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
rh=shiftdim(reshape(shiftdim(reshape(rh,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
season=shiftdim(reshape(shiftdim(reshape(season,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
day=reshape(shiftdim(day,2),365,s(1)*s(2))';
nyrs=size(rh,3);

lf=NaN*size(day);

f2=find(lat>30);f1=find(lat<-30);f=find(lat>=-10 & lat<=10);
f4=find(lat<=30 & lat>10);f5=find(lat>=-30 & lat<-10);

 %46N: Canadian standard, latitude >= 30N   (Van Wagner 1987)
  ell01 =[6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 8, 7, 6];
  %For 30 > latitude >= 10
  ell02 =[7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1,8.6, 8.1, 7.8];
  
  %#20S: For -10 > latitude >= -30  
  ell03 =[10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 9.4, 9.9, 10.2];
  %#40S: For -30 > latitude
  ell04 =[11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 10, 11.2, 11.8];
  %#For latitude near the equator, we simple use a factor of 9 for all months
 
for mi=1:12
    d1=datenum(2013,mi,1)-datenum(2013,1,0);
    d2=datenum(2013,mi+1,1)-datenum(2013,1,0)-1;
    lf(f2,d1:d2)=ell01(mi);
    lf(f1,d1:d2)=ell04(mi);
    lf(f4,d1:d2)=ell02(mi);
    lf(f5,d1:d2)=ell03(mi);
    lf(f,d1:d2)=9;
end

%constrain low end of temperature
tmx=max(tmx,-1.1);

%Eq. 16 - The log drying rate
for ii=1:nyrs
rk(:,:,ii)=squeeze(1.894*(tmx(:,:,ii)+1.1).*(100-rh(:,:,ii))).*lf*1e-4;
end

%  #Eq. 11 - Net rain amount
  rw = 0.92 * pt - 1.27;
  rw=max(0,rw);
  % initialize DMC with 15 
DMCY=10*ones(size(pt,1),1);lastseasony=DMCY*0;
for yr=1:size(pt,3);  % loop over yr
for ii=1:size(pt,2);  % loop over days
% Alteration to Eq. 12 to calculate more accurately
  wmi = 20 + 280./exp(0.023.*DMCY);
% #Eqs. 13a, 13b, 13c
  f=find(DMCY<=33);f2=find(DMCY>33 & DMCY<=65);f3=find(DMCY>65);
  B(f)=100./(0.5 + 0.3 .* DMCY(f));
  B(f2)=14 - 1.3 .* log(DMCY(f2));
  B(f3)=6.2 .* log(DMCY(f3)) - 17.2;
  B=reshape(B,size(DMCY));
% #Eq. 14 - Moisture content after rain
  wmr = wmi + 1000 * rw(:,ii,yr)./(48.77 + B .* rw(:,ii,yr)); 
% #Alteration to Eq. 15 to calculate more accurately
  pr0 = 43.43 .* (5.6348 - log(wmr - 20));
%  #Constrain P
pr=zeros(size(pr0));
  f=find(pt(:,ii,yr)<=1.5);pr(f)=DMCY(f);
  f=find(pt(:,ii,yr)>1.5);pr(f)=pr0(f);

  pr=max(0,pr);
  pr=reshape(pr,size(pr0));
% #Calculate final P (DMC)
  dmc1 = pr + rk(:,ii,yr);
  dmc1=max(0,dmc1);
  dmc(:,ii,yr)=dmc1;
  f=find(season(:,ii,yr)==0);
  if ~isempty(f)
  dmc(f,ii,yr)=0;
  DMCY(f)=0;
  end
  f=find(lastseasony==0 & season(:,ii,yr)==1);DMCY(f)=10; % standard season start DC no-overwintering
  lastseasony=season(:,ii,yr);
  DMCY=dmc(:,ii,yr);
end
end


