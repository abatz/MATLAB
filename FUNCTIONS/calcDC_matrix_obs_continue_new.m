function [dc]=calcDC(pt,tmx,lat, season);
%tmx=tmx-273.15;
%pt=pt*25.4;
s=size(pt);
if s(1)==1|s(2)==1

k1=repmat((1:365)',[1 size(tmx,1)]);
k2=repmat(lat',[365 1]);

day=real(day_length(k1,k2));
else
day=day_length(repmat((1:365)',[1 size(tmx,1) size(tmx,2)]),shiftdim(repmat(lat,[1 size(tmx,2) 365]),2));
end
pt=shiftdim(reshape(shiftdim(reshape(pt,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
tmx=shiftdim(reshape(shiftdim(reshape(tmx,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);
season=shiftdim(reshape(shiftdim(reshape(season,s(1),s(2),365,s(3)/365),2),365,s(3)/365,s(1)*s(2)),2);

day=reshape(shiftdim(day,2),365,s(1)*s(2))';

lf=NaN*size(day);

f1=find(lat>20);f2=find(lat<-20);f=find(lat>=-20 & lat<=20);
lf(f,:)=1.4;
  fl01 =[-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4, -1.6, -1.6];
  %20S: South of 20 degrees S
  fl02 =[6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 0.9, 3.8, 5.8];

for mi=1:12
    d1=datenum(2013,mi,1)-datenum(2013,1,0);
    d2=datenum(2013,mi+1,1)-datenum(2013,1,0)-1;
    lf(f1,d1:d2)=fl01(mi);
    lf(f2,d1:d2)=fl02(mi);
    lf(f,d1:d2)=1.4;
end
%constrain temp to -2.8
 tmx =max(tmx,-2.8);
  
%#Eq. 22 - Potential Evapotranspiration
  pe = (0.36 * (tmx + 2.8) + repmat(lf,[1 1 size(tmx,3)])) / 2;

% effective rainfall is rainfall in mm times 0.83 minus 1.27
er=pt*0.83-1.27;
er=max(0,er);

  
% initialize DC with 15 
DCY=15*ones(size(pt,1),1);lastseasony=DCY*0;lastfallDC=DCY*0;winterprecip=DCY*0;
for yr=1:size(pt,3);  % loop over yr
for ii=1:size(pt,2);  % loop over days
   smi=800*exp(-1*DCY/400);
   dr0=DCY-400*log(1+3.937*er(:,ii,yr)./smi);
   dr0=max(dr0,0);
   dr=dr0;
   f=find(pt(:,ii,yr)<=2.8);
   dr(f)=DCY(f);
   dr=reshape(dr,size(dr0));
   ddc=dr+pe(:,ii,yr);
   ddc=max(ddc,0);
   dc(:,ii,yr)=ddc;
   DCY=dc(:,ii,yr);
   f=find(season(:,ii,yr)==0);
   dc(f,ii,yr)=0;DCY(f)=0;
   f=find(lastseasony==0 & season(:,ii,yr)==1);
%%   DCY(f)=15; % standard season start DC no-overwintering
   DCY(f)=overwinterDC(lastfallDC(f),winterprecip(f));
   f2=find(lastseasony==1 & season(:,ii,yr)==0); % it is the last fall day, set lastfallDC and reset precipitation counter
   lastfallDC(f2)=DCY(f2);winterprecip(f2)=0;
   f3=find(lastseasony==0 & season(:,ii,yr)==0); % it is still winter, accumulate precipitation
   winterprecip(f2)=winterprecip(f2)+pt(f2,ii,yr);
   lastseasony=season(:,ii,yr);
   
end
end


