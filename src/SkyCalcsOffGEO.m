% Geographic calculations -- requires Mapping toolbox

close all;

% Constants
Re = 6371000;
Rgeo = 42164000;
Ageo = 35786000;
Mu = 3.986e14;
Vgeo = sqrt(Mu/Rgeo);
E = wgs84Ellipsoid('meter');

% Locations

scopelat = dms2degrees([19 49 35]); % Keck, HI
% scopelon = dms2degrees([-155 28 27]);
scopelon = 70;
scopealt = 4145;
scopename = 'Keck';

% scopelat = dms2degrees([31 41 18]); % MMT, AZ
% % scopelon = dms2degrees([-110 53 6]);
% scopelon = 70;
% scopealt = 2616;
% scopename = 'MMT';

% scopelat = -30.24073; % Giant Magellan, Chile
% % scopelon = -70.73659;
% scopelon = 70;
% scopealt = 2722;
% scopename = 'GMT';

[scopex,scopey,scopez] = geodetic2ecef(E,scopelat,scopelon,scopealt);

lgslatvec = -90:90;
lgslonvec = -180:180;

[lgslat,lgslon] = ndgrid(lgslatvec,lgslonvec);
[lgsx,lgsy,lgsz] = geodetic2ecef(E,lgslat,lgslon,Ageo);

dvs = abs((pi/2)*deg2rad(lgslat)*Vgeo);

dx = lgsx-scopex;
dy = lgsy-scopey;
dz = lgsz-scopez;

dec = rad2deg(atan2(dz,sqrt(dx.^2+dy.^2))); % declination
rtas = rad2deg(atan2(dy,dx)); % Right ascension

[azs,els,ras] = geodetic2aer(lgslat,lgslon,Ageo,scopelat,scopelon,scopealt,E);

figureMap = figure;

% Hubble Deep Field (north), HDF South, HU(X)DF/Chandra South
deeplons = [189.2058,338.2343,53.1625];
deeplats = [62.2161,-60.5507,-27.7914];

axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
title(sprintf('Delta-V cost to observe LUVOIR targets from %s',scopename))
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
% [Cel,hel] = contourm(lgslat,lgslon,els);

% [Cel,hel] = contourm(dec,rtas,els);
% clabelm(Cel,hel);

[Cel,hel] = contourm(dec,rtas,els,[10 10], 'linewidth', 2,'LineColor',[1 0 0]);

[Cdv,hdv] = contourm(dec,rtas,dvs, 'linewidth', 2);

tv = clabelm(Cdv,hdv);
tv.set('FontSize',14);
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,sprintf('SkyMap_offGEO %s.png',scopename))