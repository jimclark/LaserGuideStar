% Power needs
% 4x 6U panels and "2x2U" panels from GOMSpace
maxPwrGen = 1.15*(4*16+2*5); % 85 W

maxPwrThrust = 0.3*75 + 0.7*30; % Power cycle of the thruster
maxPwrLaser = 30;
maxPwrDraw = maxPwrLaser+ maxPwrThrust;

typPwrGen = 0.7*maxPwrGen;

pwrDef = maxPwrDraw-typPwrGen;

pwr_cap = 77*2; % Whr

pwrDur = pwr_cap/pwrDef;

pwrAng = acosd(maxPwrDraw/maxPwrGen);

% ADCS
% Sun torque

torque_solar = (solarConst/c)*(4*.2*.3+.2*.2)*cosd(45)*1.5*(0.15)*cosd(45); % Solar torque from the "flower" panels.

max_total_ang_mom = torque_solar*max(obs_dur)*daysec;
torque_thruster = speed_factor*sc_max_thrust_nom*sind(10)*0.15; % max torque from the thruster
desat_factor = torque_thruster/torque_solar;