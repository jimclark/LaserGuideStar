Tnd = 365.25*24*60*60/(2*pi);
muSE = 3.036e-6;
x1 = -muSE;
x2 = 1-muSE;

dt = (3600/Tnd);  % 1 hour steps.

yScopeInit = [1.00717285919175; 0; 0; 0; 0.0163636; 0];

tspan = 0:(3600/Tnd):3.2; % ~6 mos, time for scope to stay on-target
% tspan = 0:(60/Tnd):(1*24*3600/Tnd); % 1 day
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
xviScope = xvScope.*cos(tScope) - yvScope.*sin(tScope) - ypiScope;
yviScope = yvScope.*cos(tScope) + xvScope.*sin(tScope) + xpiScope;

xaScope = xpScope + 2.*yvScope -(1-muSE).*(xpScope-x1)./((xpScope-x1).^2 + ypScope.^2 + zpScope.^2).^(3./2) -muSE.*(xpScope-x2)./((xpScope-x2).^2 + ypScope.^2 + zpScope.^2).^(3./2);
yaScope = ypScope - 2.*xvScope -(1-muSE).*ypScope./((xpScope-x1).^2 + ypScope.^2 + zpScope.^2).^(3./2) -muSE.*ypScope./((xpScope-x2).^2 + ypScope.^2 + zpScope.^2).^(3./2);
zaScope = -(1-muSE).*zpScope./((xpScope-x1).^2 + ypScope.^2 + zpScope.^2).^(3./2) -muSE.*zpScope./((xpScope-x2).^2 + ypScope.^2 + zpScope.^2).^(3./2);

xaicScope = xaScope - xpScope - 2.*yvScope;
yaicScope = yaScope - ypScope + 2.*xvScope;

xaiScope = xaicScope.*cos(tScope) - yaicScope.*sin(tScope);
yaiScope = yaicScope.*cos(tScope) + xaicScope.*sin(tScope);

zaiScope = zaScope;

% plot(xpScope,ypScope,xpiScope,ypiScope)
%plot(xvScope,yvScope,xviScope,yviScope)
%plot(xaScope,yaScope,xaiScope,yaiScope)
% daspect([1 1 1])
% figurePVA = figure;
% plot(tScope,xpiScope,tScope,xviScope,tScope,xaiScope)
% xlabel('Time (radians)')
% ylabel('$x, \dot{x}, \ddot{x}$ (AU, AU/rad, AU/rad$^2$)','interpreter','latex')
% title('X-position, velocity, acceleration of L2 halo orbit (inertial frame)')
% saveas(figurePVA,'PVA.png')

% clc

idxStart = 1; % Start at beginning.
% idxStart = 300; % Trying to get to max turn...
% idxStart = find(ypScope==max(ypScope)); % First corner.
% idxStart = find(xpScope==max(xpScope)); % D loop.

tStart = tScope(idxStart);

goalAzI = 0;
goalAzRi = goalAzI-rad2deg(tStart);
goalEl = 0;

goalRange = range_LGS/AU;

goalXinit = xpScope(idxStart)+goalRange*cosd(goalEl)*cosd(goalAzRi);  % Exactly 100,000 km away, in the proper direction
goalYinit = ypScope(idxStart)+goalRange*cosd(goalEl)*sind(goalAzRi);
goalZinit = zpScope(idxStart)+goalRange*sind(goalEl);
%goalYinit = goalYinit + 100/AU; % Adding 100 m of position error

velXinit = xvScope(idxStart) + goalYinit - ypScope(idxStart);
velYinit = yvScope(idxStart) - goalXinit + xpScope(idxStart);
% velYinit = velYinit + 0.0001*(Tnd/AU); % Adding 0.1 mm/s of velocity error

yLGSInit = [xpScope(idxStart); ypScope(idxStart); zpScope(idxStart); ... 
    xvScope(idxStart); yvScope(idxStart); zvScope(idxStart); ...
    goalXinit; goalYinit; goalZinit; ... 
    velXinit; velYinit; zvScope(idxStart); ...
    goalAzI; goalEl; goalRange; ...
    0; 0; 0]; %

% tspan2 = tStart:(60/Tnd):(tStart+(3600*24*15/Tnd)); % 15 days, one-minute intervals.
tspan2 = tStart:(60/Tnd):(tStart+(3600*24/Tnd)); % one day, one-minute intervals.
% tspan2 = tStart:(60/Tnd):(tStart+(3600*3*24/Tnd)); % three days, one-minute intervals.

[tLGS,yLGSMat] = ode45(@cr3bpsepropCLazel3,tspan2,yLGSInit);

rpts = zeros(size(tspan2,2),3);
rptsnorm = zeros(size(tspan2,2),1);

for i = 1:size(tspan2,2)
    rpt = cr3bpsepropCLazel3rpt(tspan2(i),yLGSMat(i,:));
    rpts(i,:) = rpt;
    rptsnorm(i) = norm(rpt);
end

rptsI = rpts;
rptsI(:,1) = rpts(:,1).*cos(tLGS) - rpts(:,2).*sin(tLGS);
rptsI(:,2) = rpts(:,2).*cos(tLGS) + rpts(:,1).*sin(tLGS);

xpScope2 = yLGSMat(:,1);
ypScope2 = yLGSMat(:,2);
zpScope2 = yLGSMat(:,3);
xvScope2 = yLGSMat(:,4);
yvScope2 = yLGSMat(:,5);
zvScope2 = yLGSMat(:,6);

xpiScope2 = xpScope2.*cos(tLGS) - ypScope2.*sin(tLGS); % inertial reference frame
ypiScope2 = ypScope2.*cos(tLGS) + xpScope2.*sin(tLGS);
xviScope2 = xvScope2.*cos(tLGS) - yvScope2.*sin(tLGS);
yviScope2 = yvScope2.*cos(tLGS) + xvScope2.*sin(tLGS);

xpLGS = yLGSMat(:,7);
ypLGS = yLGSMat(:,8);
zpLGS = yLGSMat(:,9);
xvLGS = yLGSMat(:,10);
yvLGS = yLGSMat(:,11);
zvLGS = yLGSMat(:,12);

xpiLGS = xpLGS.*cos(tLGS) - ypLGS.*sin(tLGS); % inertial reference frame
ypiLGS = ypLGS.*cos(tLGS) + xpLGS.*sin(tLGS);
xviLGS = xvLGS.*cos(tLGS) - yvLGS.*sin(tLGS); % inertial reference frame
yviLGS = yvLGS.*cos(tLGS) + xvLGS.*sin(tLGS);

dxp = xpLGS-xpScope2;
dyp = ypLGS-ypScope2;
dzp = zpLGS - zpScope2;

range = sqrt(dxp.^2 + dyp.^2 + dzp.^2);

dxpi = xpiLGS - xpiScope2;
dypi = ypiLGS - ypiScope2;

angleAz = atan2(dypi,dxpi);
angleAzR = atan2(dyp,dxp);
dangAz = angleAz(2:end) - angleAz(1:end-1);
angleEl = atan2(dzp,sqrt(dxpi.^2+dypi.^2));

dyvi = yviLGS - yviScope2;

%%

close all

xshad = [1-muSE, 1.014, 1.014, 1-muSE];
yshad = [Re/AU, 1.0829e-04, -1.0829e-04, -Re/AU];
figureOVR = figure;
hold on;
patch(xshad,yshad,[0.4 0.4 0.6]);
plot(xpScope,ypScope,xpLGS,ypLGS)
scatter(1-muSE,0,'b*')
xlim([0.995 1.015])
daspect([1 1 1]);
hold off;
title('Telescope and LGS at L2 (AU, rotating frame)')
legend('Earth''s penumbra','Scope','Earth','LGS','Location','northwest')
saveas(figureOVR,'LGS-Scope-position-rotating.png');

figureOVI = figure;
hold on;
plot(((xpiLGS-xpiScope2)*AU/1000e3),(ypiLGS-ypiScope2)*AU, 'linewidth', 2)
% scatter(0,0,'b*', 'linewidth', 2)
plot([1.1*range_LGS/1e6,range_LGS/1e6,range_LGS/1e6,1.1*range_LGS/1e6],[1.1*iwa_box_rad,iwa_box_rad,-iwa_box_rad,-1.1*iwa_box_rad], 'linewidth', 2)
plot([range_LGS/1e6,1.1*range_LGS/1e6],[0, 0],':','linewidth',2)
% scatter(1-muSE,0,'b*')
% daspect([1 1 1]);
% hold off;
xlim([4.31e1,4.35e1])
ylim([-1,1])
title('LGS relative position to Scope (inertial frame)')
xlabel('Distance from telescope (1000''s of km)')
ylabel('Distance from line of sight (m)')
% legend('LGS','Scope','Goal','Location','northwest')
legend('LGS','Requirement (deep IWA, 0.25 lam/D)','Goal')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureOVI,'LGS-Scope-position-inertial.png');

figureAng = figure;
plot(tspan2.*Tnd./(24*60*60),angleAz*1e9,tspan2.*Tnd./(24*60*60),(iwa_box_rad/range_LGS)*1e9*ones(size(tspan2)))
% daspect([1 1 1]);
title('Scope-LGS vector azimuth (ecliptic plane, inertial frame)')
xlabel('Mission time [days]')
ylabel('Scope-LGS angle, nrad')
% saveas(figureAng,'LGS-Scope-Ecl-Ang-10k-CL-old.png');
saveas(figureAng,'LGS-Scope-Ecl-Ang-10k-CL-new.png');

figureTHR = figure;
plot(tspan2*Tnd/(3600*24),rptsnorm*(AU/Tnd^2)*24*1000)
title('Thrust over time (24 kg) for drift compensation')
xlabel('Mission time (days)')
ylabel('Total thrust (mN)')

figureTHR2I = figure;
hold on
plot(rptsI(:,1)*(AU/Tnd^2)*24*1000,rptsI(:,2)*(AU/Tnd^2)*24*1000)
quiver(0,0,0.5*cos(goalAzI),0.5*sin(goalAzI))
hold off
xlim([-1 1])
ylim([-1 1])
daspect([1 1 1])
title('Thrust vector (inertial space, 24 kg) for drift compensation')
legend('Thrust vector','Goal line of sight')
xlabel('X-thrust (mN)')
ylabel('Y-thrust (mN)')