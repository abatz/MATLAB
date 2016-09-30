function [snow]=runsnow_Dai(temp,precip);
%send in temperature in degrees C
%adopted from Dai et al., 2008 Table A1 (Land ANN): http://www.cgd.ucar.edu/cas/adai/papers/Dai_GRL08_Snow.vs.T.pdf
a=-48.2292;
b=0.7205;
c=1.1662;
d=1.0223;
snowf=a*(tanh(b*(temp-c))-d);
f=find(temp<-2);snowf(f)=100;
f=find(temp>6.5);snowf(f)=0;
snow=single(snowf./100.*precip);