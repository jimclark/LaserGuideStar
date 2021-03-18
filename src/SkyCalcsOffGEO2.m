% Geographic calculations -- requires Mapping toolbox and Phased Array?!?!
% (for rotx/roty/rotz)

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
scopelon = dms2degrees([-155 28 27]);
scopealt = 4145;
scopename = 'Keck';

% scopelat = dms2degrees([31 41 18]); % MMT, AZ
% scopelon = dms2degrees([-110 53 6]);
% scopealt = 2616;
% scopename = 'MMT';

% scopelat = -30.24073; % Giant Magellan, Chile
% scopelon = -70.73659;
% scopealt = 2722;
% scopename = 'GMT';

% daydeg = 0:360;
daydeg = 154.5:0.0001:157; % Zoom in on star 173

scopelons = scopelon + daydeg;
scopelons = scopelons - 360*(scopelons>180);

scopelats = scopelat*ones(size(daydeg));
scopealts = scopealt*ones(size(daydeg));

[scopex,scopey,scopez] = geodetic2ecef(E,scopelats,scopelons,scopealts);

lgsV0 = 68.5; % true anomaly at epoch (i.e. start of day)
lgsVs = lgsV0 + daydeg;
lgsVs = lgsVs - 360*(lgsVs>360);

lgsIP = Rgeo*[cosd(lgsVs);...
    sind(lgsVs);...
    zeros(size(daydeg))]; % Column vectors representing in-plane coordinates of LGS.

lgsinc = 15; % inclination in degrees
lgsRAAN = 32.1; % lgs RAAN, degrees

R1 = rotx(lgsinc);
R2 = rotz(lgsRAAN);

lgsvecs = R2*(R1*lgsIP);

lgsx = lgsvecs(1,:);
lgsy = lgsvecs(2,:);
lgsz = lgsvecs(3,:);

% plot3(lgsx,lgsy,lgsz)

% dvs = abs(deg2rad(lgslat)*Vgeo);

dx = lgsx-scopex;
dy = lgsy-scopey;
dz = lgsz-scopez;

decs = rad2deg(atan2(dz,sqrt(dx.^2+dy.^2))); % declination
rtas = rad2deg(atan2(dy,dx)); % Right ascension

[azs,els,ras] = geodetic2aer(lgslat,lgslon,Ageo,scopelat,scopelon,scopealt,E);

[decs_mat,starlats_mat] = ndgrid(decs,starlats);
[rtas_mat,starlons_mat] = ndgrid(rtas,starlons);
[seps,~] = distance(decs_mat,rtas_mat,starlats_mat,starlons_mat);
% [closest_to_each_tgt,idx_to_each_tgt] = min(seps,[],2);

% disp(min(min(seps)))

%%
figureApproach = figure;
plot(((daydeg-daydeg(1))*24*60*60/360)*(366.25/365.25),seps(:,173), 'linewidth', 2)
title('Angular sep. of GEO LGS from target star')
xlabel('Time in encounter (sec)')
ylabel('Angle separation (deg)')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureApproach,'GroundScope_OffGEO_window.png')

%%

figureApproachZoom = figure;
hold on
plot(((daydeg-daydeg(1))*24*60*60/360)*(366.25/365.25),seps(:,173)*3600, 'linewidth', 2)
plot([272 282],[60 60],'-.', 'linewidth', 2)
plot([272 282],[25 25],'--', 'linewidth', 2)
hold off
legend('GEO LGS separation','Keck max distance to ground LGS','Keck max distance to NGS')
ylim([0 70])
title('Angular sep. of GEO LGS from target star')
xlabel('Time in encounter (sec)')
ylabel('Angle separation (arcsec)')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureApproachZoom,'GroundScope_OffGEO_window_1arcmin.png')

%%


daydeg = 0:0.1:360;
% daydeg = 154.5:0.0001:157; % Zoom in on star 173

scopelons = scopelon + daydeg;
scopelons = scopelons - 360*(scopelons>180);

scopelats = scopelat*ones(size(daydeg));
scopealts = scopealt*ones(size(daydeg));

[scopex,scopey,scopez] = geodetic2ecef(E,scopelats,scopelons,scopealts);

lgsV0 = 68.5; % true anomaly at epoch (i.e. start of day)
lgsVs = lgsV0 + daydeg;
lgsVs = lgsVs - 360*(lgsVs>360);

lgsIP = Rgeo*[cosd(lgsVs);...
    sind(lgsVs);...
    zeros(size(daydeg))]; % Column vectors representing in-plane coordinates of LGS.

lgsinc = 15; % inclination in degrees
lgsRAAN = 32.1; % lgs RAAN, degrees

R1 = rotx(lgsinc);
R2 = rotz(lgsRAAN);

lgsvecs = R2*(R1*lgsIP);

lgsx = lgsvecs(1,:);
lgsy = lgsvecs(2,:);
lgsz = lgsvecs(3,:);

% plot3(lgsx,lgsy,lgsz)

% dvs = abs(deg2rad(lgslat)*Vgeo);

dx = lgsx-scopex;
dy = lgsy-scopey;
dz = lgsz-scopez;

decs = rad2deg(atan2(dz,sqrt(dx.^2+dy.^2))); % declination
rtas = rad2deg(atan2(dy,dx)); % Right ascension

[decs_mat,starlats_mat] = ndgrid(decs,starlats);
[rtas_mat,starlons_mat] = ndgrid(rtas,starlons);
[seps,sep_azs] = distance(decs_mat,rtas_mat,starlats_mat,starlons_mat);
[closest_to_each_tgt,idx_to_each_tgt] = min(seps,[],2);

close_candidates = starids(unique(idx_to_each_tgt(closest_to_each_tgt<0.5)));

figureMap = figure;

axesm('MapProjection','robinson','Grid','on','GLineWidth',2,'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime')
title(sprintf('Line of sight from %s through inclined GEO LGS',scopename))
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
p4 = plotm(decs,rtas,'linewidth',2);
legend([p1 p3 p2 p4],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields','LGS orbit trace'})

set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,sprintf('SkyMap_offGEO_track %s.png',scopename))