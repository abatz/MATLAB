function climograph(temp,ppt,tHIGH,tLOW,pHIGH);
% function climograph(temp,ppt,tHIGH,tLOW,pHIGH);
%Creates climograph for standard monthly temperature and precipitation
% default is for temperature in degrees C and precip in mm
% optimal inputs:
%  	tHIGH: pre-defined upper bounds for temperature
%  	tLOW:  pre-defined lower bounds for temperature
%  	pHIGH: pre-defined upper bounds for precipitation

cla;g=bar(1:1:12,ppt);
set(g(1),'FaceColor','g');
ylabel('Precipitation (mm)','fontsize',16);
h1=gca;hold on;
set(gca,'fontsize',16);
axis tight
if nargin==5
g=axis;
axis([g(1) g(2) g(3) pHIGH]);
end
set(gca,'xtick',0,'xticklabel','');
set(gca,'YAxisLocation','right');
h2 = axes('Position',get(h1,'Position'));
plot(1:1:12,temp,'r','LineWidth',3)
ylabel('Temperature (\circC)','fontsize',16);
set(h2,'YAxisLocation','left','Color','none','XTickLabel',[])
set(h2,'XLim',get(h1,'XLim'),'Layer','top')
set(gca,'xtick',1:1:12,'xticklabel','J|F|M|A|M|J|J|A|S|O|N|D');
set(gca,'fontsize',16);%axis tight
if nargin==5
g=axis;axis([g(1) g(2) tLOW tHIGH])
end
