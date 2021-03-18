%% Take a slice at constant elevation, sweep over azimuth and time of year.

AU = 1.496e11;
Re = 6371000;
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

% idxStart = 1; % Start at beginning.
% idxStart = 300; % Trying to get to max turn...
% idxStart = find(ypScope==max(ypScope)); % First corner.
% idxStart = find(xpScope==max(xpScope)); % D loop.

% disp(tStart*365.25/(2*pi))

azvec = 0:0.5:360;
tidxvec = 1:24:(1+24*180);  % Index of time points

% azvec = 92:0.001:96; % Super high resolution!
% tidxvec = 1:25;  % Once per hour for a day

[tidx,azimuths] = ndgrid(tidxvec,azvec);

accs = zeros(size(azimuths));
thrusts = zeros(size(azimuths));
times = zeros(size(azimuths));

for i = 1:numel(azimuths)

idxStart = tidx(i);
tStart = tScope(idxStart);
times(i) = tStart;

goalAzI = azimuths(i);
goalAzRi = goalAzI-rad2deg(tStart);
goalEl = 0;

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

dposin = norm(dposi);
daccin = norm(dacci);

dotp = dot(dposi,dacci);

angle = acosd(dotp/(dposin*daccin));

accreq = daccin*sind(angle);

accreqSI = accreq*(AU/Tnd^2);
TreqSI = accreqSI*sc_mass_opt_tot;

accs(i) = accreqSI;
thrusts(i) = TreqSI;

end

figureMap = figure;
[Cont,handle] = contour(azimuths,times*(Tnd/(24*60*60)),thrusts*1000,[0.03 0.1 0.2 0.3 0.5 0.7 1 3],'linewidth',2); % 0, 15, 30, 45, 60, 75 deg
% [C,h] = contour(azimuths,times*(Tnd/(24*60*60)),thrusts*1000); % 89 deg
clabel(Cont,handle,'FontSize',14);
title(sprintf('Thrust req. hold pointing (mN, elev = %d deg, %.2g,000 km)',goalEl,range_LGS/1e6))
xlabel('Azimuth (deg)')
ylabel('Mission elapsed time (days)')
set(gca, 'fontsize', 14,'linewidth',2)

saveas(figureMap,sprintf('Cost-of-watching-thrust-time-e%d-fixed.png',goalEl));