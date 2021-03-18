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
Tgeo = 2*pi*sqrt(Rgeo^3/Mu);

[earthx,earthy,earthz] = sphere;
earthx = E.MeanRadius*earthx;
earthy = E.MeanRadius*earthy;
earthz = E.MeanRadius*earthz;

% Locations

% scopelat = dms2degrees([19 49 35]); % Keck, HI
% scopelon = dms2degrees([-155 28 27]);
% scopealt = 4145;
% scopename = 'Keck';

scopelat = dms2degrees([0 0 0]); % Geostationary telescope testbed
scopelon = dms2degrees([-120 0 0]);
scopealt = Ageo;
scopename = 'GeoTT';

% scopelat = dms2degrees([31 41 18]); % MMT, AZ
% scopelon = dms2degrees([-110 53 6]);
% scopealt = 2616;
% scopename = 'MMT';

% scopelat = -30.24073; % Giant Magellan, Chile
% scopelon = -70.73659;
% scopealt = 2722;
% scopename = 'GMT';

simtime = 0:10:14*24*60*60; % 6 days, every 10 seconds
% simtime = 2.2*daysec:10:2.4*daysec; % 10 seconds, from 2.2 to 2.4 days after epoch
% simtime = 53*3600:10:55.5*3600;

scopelons = scopelon + 360*(simtime/Tgeo);

scopelats = scopelat*ones(size(simtime));
scopealts = scopealt*ones(size(simtime));

[scopex,scopey,scopez] = geodetic2ecef(E,scopelats,scopelons,scopealts);


% LGS orbit, minimum sidereal case

% lgs_ap = rSIDmin*1000;
% lgs_pe = rLEO*1000;

% sidereal, periapsis at GEO

% lgs_ap = rSIDopt*1000;
% lgs_pe = Rgeo;

% the above sidereal orbits are attempting to match the motion of Earth *at
% the Equator*, let's try something a little slower

lgs_ap = rSIDopt*1000;
lgs_pe = Rgeo/2;

lgs_sma = (lgs_ap+lgs_pe)/2;
lgs_mam = sqrt(Mu/lgs_sma^3); % mean angular motion

lgs_ecc = (lgs_ap-lgs_pe)/(lgs_ap+lgs_pe);

lgsV0 = deg2rad(68.5); % true anomaly at epoch (i.e. start of simulation)
lgsE0 = atan2(sqrt(1-lgs_ecc^2)*sin(lgsV0),lgs_ecc+cos(lgsV0));
lgsM0 = lgsE0-lgs_ecc*sin(lgsE0);

lgsMs = lgsM0+lgs_mam*simtime;
lgsEs = ecc_from_mean(lgsMs,lgs_ecc);
lgsVs = atan2(sqrt(1-lgs_ecc^2)*sin(lgsEs),cos(lgsEs)-lgs_ecc);

lgsRs = lgs_sma*(1-lgs_ecc^2)./(1+lgs_ecc*cos(lgsVs));

lgsIP = [lgsRs.*cos(lgsVs);...
    lgsRs.*sin(lgsVs);...
    zeros(size(simtime))]; % Column vectors representing in-plane coordinates of LGS.

lgsAPE = 90; % argument of periapsis, degrees
lgsinc = 60; % inclination in degrees
lgsRAAN = 32.1; % lgs RAAN, degrees

R0 = rotz(lgsAPE);
R1 = rotx(lgsinc);
R2 = rotz(lgsRAAN);

lgsvecs = R2*(R1*(R0*lgsIP));

lgsx = lgsvecs(1,:);
lgsy = lgsvecs(2,:);
lgsz = lgsvecs(3,:);

dx = lgsx-scopex;
dy = lgsy-scopey;
dz = lgsz-scopez;

dists = sqrt(dx.^2+dy.^2+dz.^2);

decs = rad2deg(atan2(dz,sqrt(dx.^2+dy.^2))); % declination
rtas = rad2deg(atan2(dy,dx)); % Right ascension

ddecs = diff(decs)./diff(simtime); % degrees-per-second difference from one moment to the next
drtas = diff(rtas)./diff(simtime);

driftrate = sqrt(ddecs.^2+drtas.^2); % this is only valid near the equator, TODO improve

%%
figureDrift = figure;
plot(simtime(1:end-1)/daysec,driftrate*3600*1000,'linewidth',2);
ylim([0 35])
set(gca, 'fontsize', 14,'linewidth',2)
xlabel('Time (days)')
ylabel('Drift rate (mas/sec)')
saveas(figureDrift,sprintf('SkyMap_HEO_track_stability_total %s.png',scopename))

%%
figureStab = figure;
hold on
plot(drtas*3600*1000,ddecs*3600*1000,'linewidth',2);
plot([35,35,-35,-35,35],[35,-35,-35,35,35],'linewidth',2);
xlim([-500,500])
ylim([-500,500])
daspect([1 1 1])
set(gca, 'fontsize', 14,'linewidth',2)
hold off
saveas(figureStab,sprintf('SkyMap_HEO_track_stability_2d %s.png',scopename))

%%

obs_idx_start = find(simtime==53.3*3600);
obs_idx_end = find(simtime==55.1*3600);

obsx = [scopex(obs_idx_start) lgsx(obs_idx_start) lgsx(obs_idx_end) scopex(obs_idx_end)];
obsy = [scopey(obs_idx_start) lgsy(obs_idx_start) lgsy(obs_idx_end) scopey(obs_idx_end)];
obsz = [scopez(obs_idx_start) lgsz(obs_idx_start) lgsz(obs_idx_end) scopez(obs_idx_end)];

figureXYZ = figure;
hold on
plot3(lgsx,lgsy,lgsz,'linewidth',2)
plot3(scopex,scopey,scopez,'linewidth',2)
plot3(obsx,obsy,obsz,'linewidth',2)
surf(earthx,earthy,earthz)
legend('LGS orbit','Telescope latitude','Best observation vector','location','southeast')
daspect([1 1 1])
view(45,30)
hold off
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureXYZ,sprintf('SkyMap_HEO_orbit %s.png',scopename))

%%

rtas_obs = rtas(simtime>53.3*3600 & simtime < 55.1*3600);
decs_obs = decs(simtime>53.3*3600 & simtime < 55.1*3600);
simtime_obs = simtime(simtime>53.3*3600 & simtime < 55.1*3600);
seps_obs = sqrt((rtas_obs-mean(rtas_obs)).^2+(decs_obs-mean(decs_obs)).^2);
figureXY = figure;
hold on
thetas = 0:0.01:2*pi;
plot((rtas_obs-mean(rtas_obs))*3600-10,(decs_obs-mean(decs_obs))*3600-10,'linewidth',2)
plot(60*cos(thetas),60*sin(thetas),'-.','linewidth',2)
plot(25*cos(thetas),25*sin(thetas),'--','linewidth',2)
title('Angular sep. of HEO LGS from target star')
legend('HEO LGS separation','Keck max distance to ground LGS','Keck max distance to NGS')
xlabel('Delta-Dec (arcsec)')
ylabel('Delta-RA (arcsec)')
xlim([-80 80])
ylim([-80 80])
daspect([1 1 1])
hold off
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureXY,sprintf('SkyMap_HEO_track_zoom %s.png',scopename))

%%
figureMap = figure;

axesm('MapProjection','robinson','Grid','on','GLineWidth',2,'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime')
title(sprintf('Line of sight from %s through HEO LGS',scopename))
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2,'DisplayName','Stark 2015 targets');
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
p4 = plotm(decs,rtas,'linewidth',2);
legend([p1 p3 p2 p4],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields','LGS orbit trace'})

set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,sprintf('SkyMap_HEO_track %s.png',scopename))