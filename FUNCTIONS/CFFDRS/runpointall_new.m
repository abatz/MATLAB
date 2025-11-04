function [fwi,ffmc,isi,dc,dmc,bui]=runpointall(temp,rh,wind,ppt,lat);
season=seasontrack_all(squeeze(temp));
s2=size(ppt);
if ndims(temp)==2
p(:,1,:)=ppt;ppt=p;
p(:,1,:)=temp;temp=p;
p(:,1,:)=rh;rh=p;
p(:,1,:)=wind;wind=p;
p(:,1,:)=season;season=p;
else
 
p(:,1,:,:)=ppt;ppt=p;
t(:,1,:,:)=temp;temp=t;
r(:,1,:,:)=rh;rh=r;
w(:,1,:,:)=wind;wind=w;
s(:,1,:,:)=season;season=s;clear s
clear p t r w
s=size(ppt);
ppt=reshape(ppt,s(1),s(2),s(3)*s(4));
temp=reshape(temp,s(1),s(2),s(3)*s(4));
rh=reshape(rh,s(1),s(2),s(3)*s(4));
wind=reshape(wind,s(1),s(2),s(3)*s(4));
season=reshape(season,s(1),s(2),s(3)*s(4));
end
whos ppt temp rh lat season
ffmc=calcFFMC_matrix_obs_continue_new(ppt,temp,rh,wind,season);
dmc=calcDMC_matrix_obs_continue_new(ppt,temp,rh,lat',season);
dc=calcDC_matrix_obs_continue_new(ppt,temp,lat,season);
dc=reshape(dc,size(wind));
dmc=reshape(dmc,size(wind));
ffmc=reshape(ffmc,size(wind));
isi=calcISI_matrix(ffmc,wind);
isi=real(isi);
bui=calcBUI_matrix(dmc,dc);
fwi=calcFWI_matrix(isi,bui);
fwi=real(fwi);
dc=reshape(dc,s2);
dmc=reshape(dmc,s2);
ffmc=reshape(ffmc,s2);
fwi=reshape(fwi,s2);
isi=reshape(isi,s2);
bui=reshape(bui,s2);

