function [f,vcom,distance]=analogsearch(var1,var2,var3,lon,lat,lonp,latp,tol1,tol2,tol3,dist);

%var1, var2, var3 are three reference datasets, time in third dimension
%lon,lat are vectors of lat/lon
%lonp, latp are the indices of the reference pixel
%tol is the acceptable amount of relative error to consider as a analog
%dist is the distance in degrees between reference and target points to consider as an analog

[lon,lat]=meshgrid(lon,lat);
distance=sqrt((lon(:)-lon(latp,lonp)).^2+(lat(:)-lat(latp,lonp)).^2);

var1m=mean(var1,3);
var2m=mean(var2,3);
var3m=mean(var3,3);
var1s=std(var1,[],3);
var2s=std(var2,[],3);
var3s=std(var3,[],3);
clear var1 var2 var3

tol1=sqrt(var1s(latp,lonp)^2)*.66;
tol2=sqrt(var2s(latp,lonp)^2)*.66;
tol3=sqrt(var3s(latp,lonp)^2)*.66;

% abs is the relative difference, defined as by the difference divided by the standard deviation for the reference point
f=find(abs(var1m(:)-var1m(latp,lonp))<=tol1 & abs(var2m(:)-var2m(latp,lonp))<=tol2 & abs(var3m(:)-var3m(latp,lonp))<=tol3 & distance>dist);


% combined measure of standard deviation; here we need to normalize standard deviation since the units are all different
% the approach I used it to divide all values by the spatial standard deviation std  

v1=abs(var1s(:)-var1s(latp,lonp))/var1s(latp,lonp);
v2=abs(var2s(:)-var2s(latp,lonp))/var2s(latp,lonp);
v3=abs(var3s(:)-var3s(latp,lonp))/var3s(latp,lonp);

vcom=sqrt(v1(:).^2+v2(:).^2+v3(:).^2);


