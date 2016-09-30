function makemap(x,y,data,values,extra);
% creates nice map bounding lat/lon coordinates
% x and y are longitude and latitude
% data is the layer you are plotting
% values are the increments between bounds are ploting, you must provide
% actual values not just the range [not needed for pcolor option]
% requires installation of m_map and cbrewer
% "extra" will result in a contour, whereas entering 4 inputs results in
% gridded "pcolor" mapping

m_proj('lambert','lat',[min(y(:)) max(y(:))],'lon',[min(x(:)) max(x(:))]);


hold on;
if nargin==5
    m_contourf(x,y,double(data),values);
else
    d=diff(x(1,:));
    m_pcolor(x-d(1)/2,y-d(1)/2,double(data));
end
shading flat;
set(gca,'climmode','manual','clim',[values(1) values(length(values))]);
m_gshhs_h('color','k');
m_grid('xtick',[],'ytick',[])
set(gcf,'paperpositionmode','auto');
shading flat
