% Geographic calculations -- requires Mapping toolbox

close all;

% Constants
Re = 6371000;
Rgeo = 42164000;
Ageo = 35786000;
E = wgs84Ellipsoid('meter');

% Locations
Klat = dms2degrees([19 49 35]); % Keck, HI
Klon = dms2degrees([-155 28 27]);
Kalt = 4145;
[Kx,Ky,Kz] = geodetic2ecef(E,Klat,Klon,Kalt);

Mlat = dms2degrees([31 41 18]); % MMT, AZ
Mlon = dms2degrees([-110 53 6]);
Malt = 2616;
[Mx,My,Mz] = geodetic2ecef(E,Mlat,Mlon,Malt);

% GSlat = -30.24073; % Gemini South, Chile
% GSlon = -70.73659;
GSlat = -30.24073; % Giant Magellan, Chile
GSlon = -70.73659;
GSalt = 2722;
[GSx,GSy,GSz] = geodetic2ecef(E,GSlat,GSlon,GSalt);

Glons = -180:0.1:180;
[Gx,Gy,Gz] = geodetic2ecef(E,0,Glons,Ageo);

[Kazs,Kels,Kras] = geodetic2aer(0,Glons,Ageo,Klat,Klon,Kalt,E);
[Mazs,Mels,Mras] = geodetic2aer(0,Glons,Ageo,Mlat,Mlon,Malt,E);
[GSazs,GSels,GSras] = geodetic2aer(0,Glons,Ageo,GSlat,GSlon,GSalt,E);

Kdx = Gx-Kx;
Kdy = Gy-Ky;
Kdz = Gz-Kz;

Kang = rad2deg(atan2(Kdz,sqrt(Kdx.^2+Kdy.^2)));

Mdx = Gx-Mx;
Mdy = Gy-My;
Mdz = Gz-Mz;

Mang = rad2deg(atan2(Mdz,sqrt(Mdx.^2+Mdy.^2)));

GSdx = Gx-GSx;
GSdy = Gy-GSy;
GSdz = Gz-GSz;

GSang = rad2deg(atan2(GSdz,sqrt(GSdx.^2+GSdy.^2)));

figCombo = figure;
hold on;
plot(Glons,Kang.*(Kels>10), 'linewidth', 2)
plot(Glons,Mang.*(Mels>10), 'linewidth', 2)
plot(Glons,GSang.*(GSels>10), 'linewidth', 2)
plot([9 21.5 25 31],[0 0 0 0],'kx', 'linewidth', 2)
plot([134],[0],'ko', 'linewidth', 2)
rectangle('Position',[-120 -0.2 40 0.4],'LineWidth',2)
h=text(9,0.2,'EDRS-A','FontSize',14);
set(h,'Rotation',45);
h=text(19,-0.2,'Artemis','FontSize',14);
set(h,'Rotation',-45);
h=text(25,0.2,'Inmarsat-4A F4','FontSize',14);
set(h,'Rotation',45);
h=text(35,-0.2,'EDRS-C','FontSize',14   );
set(h,'Rotation',-45);
text(134,0.3,'EDRS-D (TBC)','FontSize',14,'HorizontalAlignment','center');
text(-100,0.5,'LCRD (TBD)','FontSize',14,'HorizontalAlignment','center');
title('Sky coverage of demo with ground telescopes and GEO target(s)')
xlabel('Longitude of GEO target [deg]')
xlim([-180 180])
ylabel('Declination of Telescope-LGS line of sight [deg]')
legend('Keck','MMT','Gemini South')
hold off;
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figCombo,'GroundScope-GEO_Sky.png')

set(figCombo,'Units','Inches');
pos = get(figCombo,'Position');
set(figCombo,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figCombo,'GroundScope-GEO_Sky.pdf','-dpdf','-r0')

figMap = figure;
hold on;
title('Lines of sight from ground telescopes through GEO LGS')
plot([0 Re*cosd(Klat) Rgeo Rgeo+0.5*(Rgeo-Re*cosd(Klat))]/1e6,[0 Re*sind(Klat) 0 -0.5*Re*sind(Klat)]/1e6, 'linewidth', 2)
plot([0 Re*cosd(Mlat) Rgeo Rgeo+0.5*(Rgeo-Re*cosd(Mlat))]/1e6,[0 Re*sind(Mlat) 0 -0.5*Re*sind(Mlat)]/1e6, 'linewidth', 2)
plot([0 Re*cosd(GSlat) Rgeo Rgeo+0.5*(Rgeo-Re*cosd(GSlat))]/1e6,[0 Re*sind(GSlat) 0 -0.5*Re*sind(GSlat)]/1e6, 'linewidth', 2)
plot(Rgeo/1e6,0,'kx', 'linewidth', 2)
plot([-1e7,7e7]/1e6,[0 0],'k:', 'linewidth', 2)
legend('Keck/TMT','MMT','GMT')
rectangle('Position',[-Re -Re 2*Re 2*Re]/1e6,'Curvature',[1 1], 'linewidth', 2)
ylim([-10 15])
xlabel('1000''s of km')
ylabel('1000''s of km')
%xlim([-10e6 50e6])
daspect([1 1 1])
hold off;
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figMap,'GroundScope-GEO_LOS.png')
set(figMap,'Units','Inches');
pos = get(figMap,'Position');
set(figMap,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(figMap,'GroundScope-GEO_LOS.pdf','-dpdf','-r0')