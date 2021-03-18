function [dydt] = cr3bpsepropCLazel3(t,y)
% non-dimensional circular restricted 3-body problem for Sun-Earth.
% y(1) = xp, y(2) = yp, y(3) = zp,  % all of Scope
% y(4) = xv, y(5) = yv, y(6) = zv
% y(7) = xp, y(8) = yp, y(9) = zp,  % all of Laser
% y(10) = xv, y(11) = yv, y(12) = zv
% y(13) = goalAz, y(14) = goalEl, y(15) = goalRange

muSE = 3.036e-6;
AU = 1.496e11;
Tnd = 365.25*24*60*60/(2*pi);
% tDays = t*Tnd/(24*60*60);
x1 = -muSE;
x2 = 1-muSE;

xpS = y(1);
ypS = y(2);
zpS = y(3);
xvS = y(4);
yvS = y(5);
zvS = y(6);
xpL = y(7);
ypL = y(8);
zpL = y(9);
xvL = y(10);
yvL = y(11);
zvL = y(12);

goalAz = y(13);
goalEl = y(14);
goalRange = y(15);

dintXI = y(16);
dintYI = y(17);
dintZI = y(18);

xpiS = xpS.*cos(t) - ypS.*sin(t); % inertial reference frame
ypiS = ypS.*cos(t) + xpS.*sin(t);

xpiL = xpL.*cos(t) - ypL.*sin(t); % inertial reference frame
ypiL = ypL.*cos(t) + xpL.*sin(t);

goalX = xpiS+goalRange*cosd(goalEl)*cosd(goalAz);
goalY = ypiS+goalRange*cosd(goalEl)*sind(goalAz);
goalZ = zpS+goalRange*sind(goalEl);

% dxpi = xpiL-xpiS;
% dypi = ypiL-ypiS;
% dzpi = zpL - zpS;

% To convert velocities, first convert to instantaneously aligned inertial
% reference frame, then remove rotation.

xviS = xvS.*cos(t) - yvS.*sin(t) - ypiS;
yviS = yvS.*cos(t) + xvS.*sin(t) + xpiS;

xviL = xvL.*cos(t) - yvL.*sin(t) - ypiL;
yviL = yvL.*cos(t) + xvL.*sin(t) + xpiL;

% dxpr = xpL - xpS;
% dypr = ypL - ypS;

dxvi = xviL-xviS;
dyvi = yviL-yviS;
dzvi = zvL - zvS;

% range = sqrt(dxpi.^2 + dypi.^2 + dzpi.^2);
% rangeH = sqrt(dxpi.^2 + dypi.^2);

% angleAz = atan2(dypi,dxpi);
% angleAzR = atan2(dypr,dxpr);
% dangleAz = (dxpi*dyvi - dypi*dxvi)/rangeH^2; % Probably not stable around poles.
% angleEl = atan2(dzpi,rangeH);
% dangleEl = (rangeH*dzvi - dzpi*(dxpi*dxvi + dypi*dyvi)/rangeH)/(range^2); % definitely not stable around poles!

% errAz = angleAz-goalAz;
% errEl = angleEl-goalEl;

% corrAz = angleAzR - (pi/2)*sign(errAz);

xaS = xpS + 2*yvS -(1-muSE)*(xpS-x1)/((xpS-x1)^2 + ypS^2 + zpS^2)^(3/2) -muSE*(xpS-x2)/((xpS-x2)^2 + ypS^2 + zpS^2)^(3/2);
yaS = ypS - 2*xvS -(1-muSE)*ypS/((xpS-x1)^2 + ypS^2 + zpS^2)^(3/2) -muSE*ypS/((xpS-x2)^2 + ypS^2 + zpS^2)^(3/2);
zaS = -(1-muSE)*zpS/((xpS-x1)^2 + ypS^2 + zpS^2)^(3/2) -muSE*zpS/((xpS-x2)^2 + ypS^2 + zpS^2)^(3/2);

xaL = xpL + 2*yvL -(1-muSE)*(xpL-x1)/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*(xpL-x2)/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);
yaL = ypL - 2*xvL -(1-muSE)*ypL/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*ypL/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);
zaL = -(1-muSE)*zpL/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*zpL/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);

% xaD = xaL - xaS;
% yaD = yaL - yaS;
zaD = zaL - zaS;

xaicS = xaS - xpS - 2*yvS;  % First, convert to an inertial reference frame that is instantaneously co-aligned with the rotating reference frame.
yaicS = yaS - ypS + 2*xvS;

xaicL = xaL - xpL - 2*yvL;
yaicL = yaL - ypL + 2*xvL;

xaiS = xaicS.*cos(t) - yaicS.*sin(t); % Then rotate to the actual angle of inertial space.
yaiS = yaicS.*cos(t) + xaicS.*sin(t);

xaiL = xaicL.*cos(t) - yaicL.*sin(t);
yaiL = yaicL.*cos(t) + xaicL.*sin(t);

xaiD = xaiL - xaiS;
yaiD = yaiL - yaiS;

losvec = [cosd(goalEl)*cosd(goalAz); cosd(goalEl)*sind(goalAz); sind(goalEl)];

dposIx = xpiL-goalX;
dposIy = ypiL-goalY;
dposIz = zpL-goalZ;

daccI = [xaiD ; yaiD ; zaD]; % Insert sensor noise here
dvelI = [dxvi ; dyvi ; dzvi];
dposI = [dposIx; dposIy; dposIz];
dintI = [dintXI ; dintYI ; dintZI];

% dvelIinline = losvec*dot(dvelI,losvec);
% dvelIperp = dvelI-dvelIinline;

dposIinline = losvec*dot(dposI,losvec);
dposIperp = dposI-dposIinline;

if mod(t,3600*24/Tnd) == 0 % Every 24 hours
    format long
    disp(dposIperp*AU);
    format short
end

kacc = 1;% + 0.01*tanh(norm(dposIperp)*AU/4); Do not change from 1!
kvel = 10;% + 9*tanh(norm(dposIperp)*AU/4); %1000; % Velocity control scaling -- ~1 seems to work best at low displacements, ~10 at moderate
kpos = 1000;% + 1000*tanh(norm(dposIperp)*AU/400); %100000; % Position control scaling, needs ~1000 to matter
kint = 100 + 1000*tanh(norm(dposIperp)*AU/4);%1000; % Integral of position control scaling.

daccIinit = -kacc*daccI - kvel*dvelI - kpos*dposI -kint*dintI;

daccIinline = losvec*(dot(daccIinit,losvec));
daccIperp = daccIinit-daccIinline;

daccIcommand = daccIperp;
% daccIcommand = daccIinit;

% Insert thruster noise here

% daccIcommand = daccIcommand.*random('Normal',1,0.01,size(daccIcommand));

xThrI = daccIcommand(1);  
yThrI = daccIcommand(2);
zThr = daccIcommand(3);

xThrR = xThrI.*cos(t) + yThrI.*sin(t);
yThrR = yThrI.*cos(t) - xThrI.*sin(t);

xaL = xaL + xThrR;
yaL = yaL + yThrR;
zaL = zaL + zThr;

dydt = [xvS; yvS; zvS; xaS; yaS; zaS; xvL; yvL; zvL; xaL; yaL; zaL; 0; 0; 0; dposIx; dposIy; dposIz];
% rpt = [xThrR;yThrR;zThr];

end