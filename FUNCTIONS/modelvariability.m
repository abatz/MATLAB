function modelvariability(hist,future);

% from Renaud Barbero
% hist and future are 2 dimensional (space x time)
% requires notBoxplot2 and rotateXLabels

 
% change the models here
MODELS={'inmcm4','CSIRO-Mk3-6-0','CanESM2','CNRM-CM5','MIROC5','GFDL-ESM2M',...
    'GFDL-ESM2G','MRI-CGCM3','HadGEM2-ES365','HadGEM2-CC365','bcc-csm1-1',...
    'MIROC-ESM','MIROC-ESM-CHEM','BNU-ESM','bcc-csm1-1-m',...
    'IPSL-CM5A-LR','IPSL-CM5B-LR'};

% labels for space here
LABELS ={'Mixed wood shield','Western cordillera',...
    'Appalachian forests','SE coastal plains',...
    'Temperate prairies','WC semi-arid prairies','SC semi-arid prairies',...
    'Cold deserts','Warm deserts',...
    'Mediterranean CA','Western Sierra Madre',...
    'Upper Gila Mountain','Everglades'};
symbols = {'o','o','o','o','v','v','v','v',...
    'square','square','square','square','*','*','*','*','v'};

%COL=cbrewer('qual','Set1',17);
fc = {'k','b','r','g','k','b','r','g',...
    'k','b','r','g','k','b','r','g','m'};
 
% axes('position',[.19 .15 .34 .7])
% future season
axes('position',[.05 .3 .4 .4])
S = 12;
h=notBoxPlot2(future',[],[],'patch',symbols,fc);
d=[h.data]; %box

set(gca,'ytick',[0:1:10],'yticklabel',[0:1:10],'fontsize',S)
ylabel('Blah blah',...
    'Interpreter','latex','fontsize',S)
set(gca,'xtick',[1:length(LABELS)],'xticklabel',LABELS,'fontsize',S)
rotateXLabels(gca, 90 )
hold on
% historical season
% plotting range + mean
for i=1:length(LABELS)
    plot([i+.0 i+.0],[min(hist(i,:)) max(hist(i,:))],'k-');
    plot([i-.3 i+.3],[mean(hist(i,:)) mean(hist(i,:))],'k-','linewidth',2);
end
 
% legend (part1)
a=legend(MODELS);
set(a,'fontsize',8,'Location','Northeast');
