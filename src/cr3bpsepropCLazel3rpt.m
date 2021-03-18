function [rpt] = cr3bpsepropCLazel3rpt(t,y)
% Extract thrust metrics 
% non-dimensional circular restricted 3-body problem for Sun-Earth.
% y(1) = xp, y(2) = yp, y(3) = zp,  % all of Scope
% y(4) = xv, y(5) = yv, y(6) = zv
% y(7) = xp, y(8) = yp, y(9) = zp,  % all of Laser
% y(10) = xv, y(11) = yv, y(12) = zv
% y(13) = goalAz, y(14) = goalEl, y(15) = goalRange

muSE = 3.036e-6;
x1 = -muSE;
x2 = 1-muSE;

xpL = y(7);
ypL = y(8);
zpL = y(9);
xvL = y(10);
yvL = y(11);

% Accelerations under CR3BP only

xaL = xpL + 2*yvL -(1-muSE)*(xpL-x1)/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*(xpL-x2)/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);
yaL = ypL - 2*xvL -(1-muSE)*ypL/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*ypL/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);
zaL = -(1-muSE)*zpL/((xpL-x1)^2 + ypL^2 + zpL^2)^(3/2) -muSE*zpL/((xpL-x2)^2 + ypL^2 + zpL^2)^(3/2);

dydt = cr3bpsepropCLazel3(t,y);

xaLt = dydt(10);
yaLt = dydt(11);
zaLt = dydt(12);

rpt = [xaLt-xaL;yaLt-yaL;zaLt-zaL];

end