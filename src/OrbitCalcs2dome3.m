%% Run a full cube of az, el, and time, average over time, divide dV cap by it, divide by 2.

AU = 1.496e11;
Re = 6371000;
Tnd = 365.25*24*60*60/(2*pi);
muSE = 3.036e-6;
aMoon = 384400e3;
aGeo = 42164000;
dt = (3600/Tnd);  % 1 hour steps.

%%


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
elvec = 0:0.5:90;
tidxvec = 1:24:numel(tScope);  % Index of time points, skipping 1 day at a time.

[elevs,azimuths,tidx] = ndgrid(elvec,azvec,tidxvec);

[elevs2d, azimuths2d] = ndgrid(elvec,azvec);



avgaccs = zeros(size(azimuths2d));

accs = zeros(size(azimuths));
thrusts = zeros(size(azimuths));
times = zeros(size(azimuths));

desrange = range_LGS/AU;

%%

total = numel(azimuths);

for i = 1:total
    
    
    if mod(i,10000) == 0
        clc
        fprintf('%.1f%%\n',100*i/total)
    end
    
    
    idxStart = tidx(i);
    tStart = tScope(idxStart);
    times(i) = tStart;
    
    goalAzI = azimuths(i);
    goalAzRi = goalAzI-rad2deg(tStart);
    goalEl = elevs(i);
    
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


%%
for i = 1:numel(elvec)
    for j = 1:numel(azvec)
        avgaccs(i,j) = mean(accs(i,j,:));
    end
end

max_bg_acc = max(max(max(accs)));
min_bg_acc = min(min(min(accs)));
max_bg_thrust = max_bg_acc*sc_mass_opt_tot;
min_bg_thrust = min_bg_acc*sc_mass_opt_tot;

avg_acc = mean(mean(avgaccs));

% Technically this should go to a later script, after prop selection

timeObs = 0.5*(dvcap./avgaccs);

%%

% Contours = [100 150 200 300];
Contours = [300 600 1200 2400 4800];
figureMap = figure;
surf(azimuths2d,elevs2d, log(timeObs/(24*60*60)),'EdgeColor','none');
colorbar('YTick',log(Contours),'YTickLabel',Contours);
colormap(jet);
caxis(log([Contours(1) Contours(length(Contours))]));
colorbar('FontSize',12,'YTick',log(Contours),'YTickLabel',Contours);
title(sprintf('Max. obs. time (days, 12U RF Ion, %.2g,000 km range)',range_LGS/1e6))
xlabel('Ecliptic longitude (deg)')
ylabel('Ecliptic latitude (deg)')
view([0 90]);
xlim([0 360]);
ylim([0 90]);

% set(findall(gcf,'-property','FontSize'),'FontSize',14);

% save('cost-of-observation-10k.mat','accs','avgaccs','elevs','azimuths','tidx')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,sprintf('Max-watching-time-%.2gk.png',range_LGS/1e6));

%%

xshad = [1-muSE, 1.014, 1.014, 1-muSE];
yshad = [Re/AU, 1.0829e-04, -1.0829e-04, -Re/AU];

tOrb = 0:360;
xMoon = (1-muSE) + (aMoon/AU)*cosd(tOrb);
yMoon = (aMoon/AU)*sind(tOrb);

xGeo = (1-muSE) + (aGeo/AU)*cosd(tOrb);
yGeo = (aGeo/AU)*sind(tOrb);

figureOVR = figure;
hold on;
plot(xpScope,ypScope,'linewidth',2)
scatter(1-muSE,0,'b*','linewidth',2)
patch(xshad,yshad,[0.4 0.4 0.6],'linewidth',2);
plot(xMoon,yMoon,'linewidth',2,'color',[0.5 0.5 0.5]);
quiver(1.003,-0.005,0.0037,0.0045,0)
text(1.0027,-0.005,'t=0','FontSize',14,'HorizontalAlignment','right')
% annotation('textarrow',[1.00717285919175 1.003],[0 -0.005],'String','t = 0 ')
% plot(xGeo,yGeo,'k--','linewidth',2);
% scatter(-muSE,0,'y*')
daspect([1 1 1]);
xlim([0.995 1.015])
xlabel('AU')
ylabel('AU')
hold off;
title('Telescope orbit at L2 (rotating frame)')
legend('L2 halo orbit','Earth','Earth''s penumbra','Moon''s orbit','Location','northwest')
set(gca,'linewidth',2,'FontSize',14)
saveas(figureOVR,'Scope-position-rotating.png')

%%

figureOVR_all = figure;
hold on;
plot(xpScope,ypScope,'linewidth',2)
scatter(1-muSE,0,'b*','linewidth',2)
scatter(-muSE,0,'r*','linewidth',2)
% scatter([1-muSE-(muSE/3)^(1/3),1-muSE+(muSE/3)^(1/3),-1-muSE-(7*muSE/12),0.5-muSE,0.5-muSE],[0,0,0,sqrt(3)/2,-sqrt(3)/2],'g*','linewidth',2)
scatter([1-muSE-(muSE/3)^(1/3),1-muSE+(muSE/3)^(1/3)],[0,0],'g*','linewidth',2)
% patch(xshad,yshad,[0.4 0.4 0.6],'linewidth',2);
plot(xMoon,yMoon,'linewidth',2,'color',[0.5 0.5 0.5]);
% plot(xGeo,yGeo,'k--','linewidth',2);
% scatter(-muSE,0,'y*')
daspect([1 1 1]);
ylim([-0.1 0.1])
xlim([-0.1 1.1])
hold off;
title('Telescope orbit at L2 (AU, rotating frame)')
legend('L2 Halo orbit','Earth','Sun','L1/L2 Lagrange points','Moon''s orbit','Location','northwest')
set(gca,'linewidth',2,'FontSize',14)
saveas(figureOVR_all,'Scope-position-rotating-all-lagrange.png')