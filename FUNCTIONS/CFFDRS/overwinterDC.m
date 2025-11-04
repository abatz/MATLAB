function [lastDC]=overwinterDC(lastDC,rw,a,b);
if nargin==2 %default values
a=0.75;
b=0.75;
end
Qf = 800*exp(-lastDC/400);
Qs = a*Qf+b*(3.94*rw);
lastDC = 400*log(800./Qs);
lastDC(lastDC<15)=15;

