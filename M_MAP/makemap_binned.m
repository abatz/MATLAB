function makemap(x,y,data,values);
% creates nice map bounding lat/lon coordinates
% x and y are longitude and latitude
% data is the layer you are plotting
% values are the increments between bounds are ploting, you must provide
% actual values not just the range
% requires installation of m_map and cbrewer
cla;
m_proj('lambert','lat',[min(y(:)) max(y(:))],'lon',[min(x(:)) max(x(:))]);

data2=NaN(size(data));
f=find(data<=values(1));
data2(f)=1;
for i=1:length(values)-1
 f=find(data>values(i) & data<=values(i+1));
 data2(f)=i+1;
end
f=find(data>values(length(values)));
data2(f)=length(values)+1;

hold on;
m_pcolor(x,y,double(data2));
colormap(cbrewer('seq','YlOrRd',length(values)));
h=colormap;
h(2:size(h,1)+1,:)=h;
h(1,:)=[1 1 1];
colormap(h)
shading flat;
%m_grid;
set(gca,'climmode','manual','clim',[1 2+(length(values))]);
%m_coast('color','k');
m_gshhs('ic','color',[.6 .6 .6],'linewidth',2);
m_grid('xtick',[],'ytick',[])
set(gcf,'paperpositionmode','auto');
%set(gca, 'Box','on');
shading flat
