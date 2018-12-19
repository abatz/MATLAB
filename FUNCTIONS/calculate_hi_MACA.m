 % list of models to run calculations over
 MODEL_NAME={'inmcm4'; 'CSIRO-Mk3-6-0';'CanESM2';'CNRM-CM5';'MIROC5';'GFDL-ESM2M';'GFDL-ESM2G';'MRI-CGCM3';'HadGEM2-ES365';'HadGEM2-CC365';'bcc-csm1-1';'MIROC-ESM';'MIROC-ESM-CHEM';'BNU-ESM';'bcc-csm1-1-m';'CCSM4';'IPSL-CM5A-LR';'IPSL-CM5A-MR';'IPSL-CM5B-LR';'NorESM1-M'};

% loop over models
for ii=[1:20]

% choose time period/scenario
time='_historical';
%time='_rcp45';
%time='_rcp85';

% load pointer to tmax and rmin
dd='/data/gcm/cmip5/24thdeg/macav2metdata/daily/';
mx=matfile([dd,'maca_v2_metdata_2var_10pat_CONUS_',char(MODEL_NAME(ii)),time,'_tasmax.mat']);
rn=matfile([dd,'maca_v2_metdata_2var_10pat_CONUS_',char(MODEL_NAME(ii)),time,'_rhsmin.mat']);

% run calculations over 33 blocks of 42-latitude points, all longitudes, all days
for i=1:33
ij=(i-1)*42+1:i*42;
% convert temperature to degrees F
temp=(mx.data(:,:,:,ij)-273.15)*1.8+32;
relh=rn.data(:,:,:,ij);
% only use data from April 1 to Oct 31
temp=temp(91:304,:,:,:);
relh=relh(91:304,:,:,:);

% calculate heat index
[heat]=heatindex(temp(:),relh(:));
heat=reshape(heat,size(temp));

% next write data to file
m=matfile(['heatindex',time,'_',char(MODEL_NAME(ii))]);
m.Properties.Writable=true;

if i==1
% if this is the first time writing to the file
m.heat=uint8(heat);
else
% otherwise append file
m.heat(:,:,:,ij)=uint8(heat);
clear heat in temp relh;
end
end
clear m
end


