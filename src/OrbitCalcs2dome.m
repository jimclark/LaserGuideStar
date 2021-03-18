% At a given time, show thrust needed vs. az-el

close all
Tnd = 365.25*24*60*60/(2*pi);
muSE = 3.036e-6;

dt = (3600/Tnd);  % 1 hour steps.

yScopeInit = [1.00717285919175; 0; 0; 0; 0.0163636; 0];

tspan = 0:(3600/Tnd):3.2; % ~6 mos, time for scope to stay on-target
%tspan = 0:0.001:1; %~2 mos, time for LGS to stay close without active station-keeping
%tspan = linspace(0,24*60*60/Tnd,10000); % 1 day

[tScope,yScopeMat] = ode45(@cr3bpse,tspan,yScopeInit);

xpScope = yScopeMat(:,1);
ypScope = yScopeMat(:,2);
zpScope = yScopeMat(:,3);
xvScope = yScopeMat(:,4);
yvScope = yScopeMat(:,5);
zvScope = yScopeMat(:,6);

xpiScope = xpScope.*cos(tScope) - ypScope.*sin(tScope); % inertial reference frame
ypiScope = ypScope.*cos(tScope) + xpScope.*sin(tScope);
xviScope = xvScope.*cos(tScope) - yvScope.*sin(tScope);
yviScope = yvScope.*cos(tScope) + xvScope.*sin(tScope);

idxStart = 1; % Start at beginning.
% idxStart = 300; % Trying to get to max turn...
% idxStart = find(ypScope==max(ypScope)); % First corner, 45 days.
% idxStart = find(xpScope==max(xpScope)); % outside of D loop, 93 days.

tStart = tScope(idxStart);
tDays = tStart*365.25/(2*pi);

azvec = 0:0.5:360;
elvec = 0:0.5:90;

[elevs,azimuths] = ndgrid(elvec,azvec);

accs = zeros(size(azimuths));
thrusts = zeros(size(azimuths));

for i = 1:numel(azimuths)

goalAzI = azimuths(i);
goalAzRi = goalAzI-rad2deg(tStart);
goalEl = elevs(i);

desrange = range_LGS/AU;

goalXinit = xpScope(idxStart)+desrange*cosd(goalEl)*cosd(goalAzRi);
goalYinit = ypScope(idxStart)+desrange*cosd(goalEl)*sind(goalAzRi);
goalZinit = zpScope(idxStart)+desrange*sind(goalEl);

velXinit = xvScope(idxStart) + goalYinit - ypScope(idxStart);
velYinit = yvScope(idxStart) - goalXinit + xpScope(idxStart);

goalXi = goalXinit*cos(tStart) - goalYinit*sin(tStart); % inertial reference frame
goalYi = goalYinit*cos(tStart) + goalXinit*sin(tStart);

radInit = sqrt(goalXi^2 + goalYi^2 + goalZinit^2);

yScope = [xpScope(idxStart); ypScope(idxStart); zpScope(idxStart); ... 
    xvScope(idxStart); yvScope(idxStart); zvScope(idxStart)]; %

yLGS = [goalXinit; goalYinit; goalZinit; ... 
    velXinit; velYinit; zvScope(idxStart)]; %

dydtLGS = cr3bpse(tStart,yLGS);
dydtScope = cr3bpse(tStart,yScope);

xvScoper = dydtScope(1);
yvScoper = dydtScope(2);
zvScoper = dydtScope(3);
xaScoper = dydtScope(4);
yaScoper = dydtScope(5);
zaScoper = dydtScope(6);

xaScopeic = xaScoper - xpScope(idxStart) - 2*yvScoper; % Coaligned inertial reference frame
yaScopeic = yaScoper - ypScope(idxStart) + 2*xvScoper;

xaScopei = xaScopeic*cos(tStart) - yaScopeic*sin(tStart); % rotate
yaScopei = yaScopeic*cos(tStart) + xaScopeic*sin(tStart);

xvLGSr = dydtLGS(1);
yvLGSr = dydtLGS(2);
zvLGSr = dydtLGS(3);
xaLGSr = dydtLGS(4);
yaLGSr = dydtLGS(5);
zaLGSr = dydtLGS(6);

xaLGSic = xaLGSr - goalXinit - 2*yvLGSr; % Coaligned inertial reference frame
yaLGSic = yaLGSr - goalYinit + 2*xvLGSr;

xaLGSi = xaLGSic*cos(tStart) - yaLGSic*sin(tStart); % rotate
yaLGSi = yaLGSic*cos(tStart) + xaLGSic*sin(tStart);

dposi = [desrange*cosd(goalEl)*cosd(goalAzI) ; desrange*cosd(goalEl)*sind(goalAzI) ; desrange*sind(goalEl)];
dacci = [xaLGSi-xaScopei ; yaLGSi-yaScopei ; zaLGSr-zaScoper ];

dacciInline = dposi*dot(dposi,dacci)/dot(dposi,dposi);
dacciCross = dacci-dacciInline;

dposin = norm(dposi);
daccin = norm(dacci);

dotp = dot(dposi,dacci);

angle = acosd(dotp/(dposin*daccin));

accreq = daccin*sind(angle);
% accreq = norm(dacciCross);

accreqSI = accreq*(AU/Tnd^2);
TreqSI = accreqSI*sc_mass_opt_tot;

accs(i) = accreqSI;
thrusts(i) = TreqSI;

end


figureMap = figure;
colormap cool
axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
[Cont,handle] = contourm(elevs,azimuths,thrusts*1000,[0.03 0.1 0.3 0.5 0.7 1 3],'linewidth',2); % 0 days, 100k km
[Cont2,handle2] = contourm(-elevs,azimuths,thrusts*1000,[0.03 0.1 0.3 0.5 0.7 1 3],'linewidth',2); % 0 days, 100k km
% [C,h] = contour(azimuths,elevs,thrusts*1000,[0.01 0.05 0.09 0.13]); % 0 days, 10k km sep
% [C,h] = contour(azimuths,elevs,thrusts*1000,[0.01 0.05 0.09 0.13]); % 45/90 days, 10k km sep
% [C,h] = contour(azimuths,elevs,thrusts*1000);
% [C,h] = contour(azimuths,elevs,thrusts*1000,[0.03 0.1 0.2 0.3 0.5 0.7 1 3]); % 90 days, D loop, 100k km sep
% clabel(C,h,'FontSize',14);
htext = clabelm(Cont,handle);
set(htext,'fontsize',14);
title(sprintf('Thrust req. hold pointing (mN, t = %d days, %.2g,000 km)',round(tDays),range_LGS/1e6))
colorbar
% xlabel('Ecliptic longitude (deg)')
% ylabel('Ecliptic latitude (deg)')
set(gca, 'fontsize', 14,'linewidth',2)

saveas(figureMap,sprintf('Cost-of-watching-thrust-t0-fixed-%.2gk.png',range_LGS/1e6));
% colormap default

% figureSphere = figure;
% xthr = (2.7e-4-thrusts).*cosd(elevs).*cosd(azimuths);  % Exactly 100,000 km away, in the proper direction
% ythr = (2.7e-4-thrusts).*cosd(elevs).*sind(azimuths);
% zthr = (2.7e-4-thrusts).*sind(elevs);
% plot3(xthr,ythr,zthr)
% daspect([1 1 1])
% 
% set(findall(gcf,'-property','FontSize'),'FontSize',14);



