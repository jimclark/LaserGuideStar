Tnd = 365.25*24*60*60/(2*pi);

yScopeInit = [1.00717285919175; 0; 0; 0; 0.0163636; 0];

tspan = 0:(3600/Tnd):3.2; % ~6 mos, time for scope to stay on-target

[tScope,yScopeMat] = ode45(@cr3bpse,tspan,yScopeInit);

xpScope = yScopeMat(:,1);
ypScope = yScopeMat(:,2);
zpScope = yScopeMat(:,3);
xvScope = yScopeMat(:,4);
yvScope = yScopeMat(:,5);
zvScope = yScopeMat(:,6);

yLGSInit = [xpScope(1); ypScope(1); zpScope(1); ... 
    xvScope(1)+1.44*(Tnd/AU); yvScope(1); zvScope(1)]; %adding pushoff velocity

[tLGS,yLGSMat] = ode45(@cr3bpse,tspan,yLGSInit);

xpLGS = yLGSMat(:,1);
ypLGS = yLGSMat(:,2);
zpLGS = yLGSMat(:,3);
xvLGS = yLGSMat(:,4);
yvLGS = yLGSMat(:,5);
zvLGS = yLGSMat(:,6);

dist = sqrt((xpLGS-xpScope).^2+(zpLGS-zpScope).^2+(ypLGS-ypScope).^2);
dvs = sqrt((xvLGS-xvScope).^2+(zvLGS-zvScope).^2+(yvLGS-yvScope).^2);

drift_thresh = find(dist*AU>range_LGS);

drift_time = Tnd*tspan(drift_thresh(1))/daysec;
drift_vel = dvs(drift_thresh(1))*AU/Tnd;

decel_time = (drift_vel*sc_mass_opt_tot/sc_max_thrust_nom)/daysec;