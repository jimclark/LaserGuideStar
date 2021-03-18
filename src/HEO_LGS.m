close all
mu = 398600;

RE = 6378;
rLEO = RE+400;
rGEO = RE+35786;
aLUN = 384400;
rLUN1 = 362600;
rLUN2 = 405400;
bLUN = sqrt(rLUN1*rLUN2);
vLEO = sqrt(mu/rLEO);
vESC = sqrt(2*mu/rLEO);
vGEO = sqrt(mu/rGEO);
tGEO = 2*pi*sqrt(rGEO^3/mu); % 1 sidereal day

rVABi1 = (1+1)*RE; % inner Van Allen Belt starts at ~1,000 km altitude (i.e. 1.17 RE from center), but it peaks at 1-3 RE altitude = 2-4 from center.
rVABi2 = (1+3)*RE;

rVABo1 = (1+4)*RE; % outer Van Allen Belt peaks 4-6 RE altitude
rVABo2 = (1+6)*RE;

vGTO1 = sqrt(2*mu*(1/rLEO-1/(rLEO+rGEO)));

vGTO2 = sqrt(2*mu*(1/rGEO-1/(rLEO+rGEO)));

aGTO = (rLEO+rGEO)/2;
bGTO = sqrt(rLEO*rGEO);

dV_leo_geo = (vGTO1-vLEO)+(vGEO-vGTO2);



vSID = 2*pi*RE/tGEO;

rSIDcirc = mu/vSID^2;


%minimum sidereal-apogee orbit, so 400 x 150,000 km orbit

rSIDmins = roots([1 rLEO -2*mu*rLEO/vSID^2]);

rSIDmin = max(rSIDmins);

tSIDmin = 2*pi*sqrt(((rSIDmin+rLEO)/2)^3/mu);

tGTO = 2*pi*sqrt(aGTO^3/mu);

vSIDminP = sqrt(2*mu*(1/rLEO-1/(rLEO+rSIDmin)));

dv_leo_sid_min = vSIDminP-vLEO;

dv_gto_sid_min = vSIDminP-vGTO1;

aSIDmin = (rLEO+rSIDmin)/2;
bSIDmin = sqrt(rLEO*rSIDmin);

max_burn_d = 2*RE/vSIDminP; % Assume spacecraft moves at a constant v at periapsis, how long does it take to move across the Earth (i.e. "time close enough to periapsis")

num_burns = (dv_gto_sid_min*1000*sc_mass_opt_tot/sc_max_thrust_nom)/max_burn_d; % Need to multiply by 1000 to get from km/s to m/s

approx_depl_duration = num_burns*(tGTO+tSIDmin)/2;

%%
% Let's calculate the deployment duration more exactly

test_sc_vel = vGTO1;
depl_duration = 0;
burn_count = 0;

while test_sc_vel < vSIDminP
    burn_count = burn_count + 1;
    burn_d = 2*RE/test_sc_vel;% Assume spacecraft moves at a constant v at periapsis, how long does it take to move across the Earth (i.e. "time close enough to periapsis")
    test_sc_vel = test_sc_vel + burn_d*sc_max_thrust_nom/(sc_mass_opt_tot*1000); % need to divide by 1000 to get from m/s to km/s
    test_sc_sma = (mu)/((2*mu/rLEO)-test_sc_vel^2);
    test_sc_pd = 2*pi*sqrt(test_sc_sma^3/mu);
    depl_duration = depl_duration+test_sc_pd;
end

% End result turns out to be 5168 burns (vs. 5300 preducted above), total
% duration 14.6 years vs. 20 predicted.  Not bad, but good to know to feed
% SPENVIS.

mass_shield = 0.001*(2*20*20+4*20*30)*0.7*2.7; % aluminum 0.7 cm (!) thick

%% What if we wanted to get there in, say, six months?
% There should be a way to solve for this directly...

test_sc_vel = vGTO1;
depl_duration_ht = 0;
burn_count_ht = 0;
thrust_factor = 30; % 30 turns out to be about right, total duration 181 days

while test_sc_vel < vSIDminP
    burn_count_ht = burn_count_ht + 1;
    burn_d = 2*RE/test_sc_vel;% Assume spacecraft moves at a constant v at periapsis, how long does it take to move across the Earth (i.e. "time close enough to periapsis")
    test_sc_vel = test_sc_vel + thrust_factor*burn_d*sc_max_thrust_nom/(sc_mass_opt_tot*1000); % need to divide by 1000 to get from m/s to km/s
    test_sc_sma = (mu)/((2*mu/rLEO)-test_sc_vel^2);
    test_sc_pd = 2*pi*sqrt(test_sc_sma^3/mu);
    depl_duration_ht = depl_duration_ht+test_sc_pd;
end


%%
rad_count = 0.05*(RE/vSIDminP)*2*num_burns; % Assume Van Allen Belt's dose is concentrated at 0.05 rad/sec over a 1 RE thickness, per: https://spacemath.gsfc.nasa.gov/Algebra1/3Page7.pdf

% Verifying radiation dosage information
vGTOvab = sqrt(mu*(2/(2*RE)-(1/aGTO))); % velocity of spacecraft in GTO at van Allen belt (i.e. 1 RE altitude, 2 RE radius)
rad_gto_yr = 0.05*2*(RE/vGTOvab)*(yrsec/tGTO); % 2 passes through the belts per orbit, times ~800 orbits per year, comes out to 100 krad per yer -- but other sources suggest more like 2.5 per year?

vSIDminvab = sqrt(mu*(2/(2*RE)-(1/aSIDmin))); % velocity of spacecraft on minimal sidereal orbit at van Allen belt (i.e. 1 RE altitude, 2 RE radius)

rad_gto_smad_noshld_yr = (3e6/yrsec)*2*(RE/vGTOvab)*(yrsec/tGTO); % 150 krad/yr, using SMAD fig 7-11 p 135.
rad_gto_smad_3mm_yr = (3e4/yrsec)*2*(RE/vGTOvab)*(yrsec/tGTO); % 1.5 krad/yr w/ 0.8 g/cm2 of aluminum everywhere (so 3 mm thick -- that's 2.5 kg to cover a whole 12U satellite)

rad_sid_smad_noshld = (3e6/yrsec)*2*(RE/vSIDminvab)*num_burns; % 800 krad!!!
rad_sid_smad_maxshld = (1e4/yrsec)*2*(RE/vSIDminvab)*num_burns; % 2 krad, even after cladding the spacecraft in 8 kg of aluminum (that's 1-cm-thick aluminum, and also 1/3rd the mass budget)

%What if we wanted to keep the periapsis at GEO?

rSIDopts = roots([1 rGEO -2*mu*rGEO/vSID^2]);

rSIDopt = max(rSIDopts); % turns out to be ~lunar distance

aSIDopt = (rGEO+rSIDopt)/2;
bSIDopt = sqrt(rSIDopt*rGEO);

tSIDopt = 2*pi*sqrt(((rSIDopt+rGEO)/2)^3/mu);

vSIDoptP = sqrt(2*mu*(1/rGEO-1/(rGEO+rSIDmin)));

% going from GEO to SID-OPT orbit
dv_geo_sid_opt = vSIDoptP-vGEO;


vSTO1 = sqrt(2*mu*(1/rLEO-1/(rLEO+rSIDopt)));

vSTO2 = sqrt(2*mu*(1/rSIDopt-1/(rLEO+rSIDopt)));

dv_gto_sid_opt = (vSTO1-vGTO1) + (vSID-vSTO2);

dv_inc = (pi/2)*vGEO*deg2rad(15);

time_inc = (dv_inc*1000*sc_mass_opt_tot/sc_max_thrust_nom);

% Make inner and outer boundaries of each belt
t = linspace(0,2*pi,100);
xi1 = rVABi1*cos(t);
xi2 = rVABi2*cos(t);
yi1 = rVABi1*sin(t);
yi2 = rVABi2*sin(t);

xo1 = rVABo1*cos(t);
xo2 = rVABo2*cos(t);
yo1 = rVABo1*sin(t);
yo2 = rVABo2*sin(t);


%%
figureOrbits = figure;
hold on
axis equal
rectangle('Position',[-RE -RE 2*RE 2*RE],'Curvature',[1 1],'FaceColor',[0 .5 .5],'LineStyle','none'); % Earth
% rectangle('Position',[-rLEO -rLEO 2*rLEO 2*rLEO],'Curvature',[1 1],'LineWidth',2);
rectangle('Position',[-rLEO -bGTO 2*aGTO 2*bGTO],'Curvature',[1 1],'LineWidth',2,'LineStyle',':'); % GTO
rectangle('Position',[-rGEO -rGEO 2*rGEO 2*rGEO],'Curvature',[1 1],'LineWidth',2,'LineStyle','--'); % GEO
title('Comparison of different LGS orbits')
xlabel('Distance from Earth center (km)')
ylabel('Distance from Earth center (km)')

plot([(aGTO-rLEO) 1.3e5],[-bGTO -1.85e5],'k:','linewidth',1.5)
plot([0 -0.7e5],[rGEO 1.35e5],'k--','linewidth',1.5)

plot([-3*RE -0.9e5 -6*RE],[0 -1e5 0],'k')

text(-1e5,1.5e5,'GEO','FontSize',14)
text(-1e5,-1e5,'Van Allen belts','FontSize',14,'HorizontalAlignment','right')
text(1e5,-2e5,'GTO','FontSize',14)

text(3e5,2.3e5,'Moon','FontSize',14)
text(1e5,0.5e5,'Sidereal orbit 1','FontSize',14)
text(1e5,1.5e5,'Sidereal orbit 2','FontSize',14)

set(gca, 'fontsize', 14,'linewidth',2)

rectangle('Position',[-rLEO -bSIDmin 2*aSIDmin 2*bSIDmin],'Curvature',[1 1],'LineWidth',2);
rectangle('Position',[-rGEO -bSIDopt 2*aSIDopt 2*bSIDopt],'Curvature',[1 1],'LineWidth',2);
rectangle('Position',[-rLUN2 -bLUN 2*aLUN 2*bLUN],'Curvature',[1 1],'LineWidth',2,'EdgeColor',[0.5 0.5 0.5]);

hp1 = patch([xi2,xi1],[yi2,yi1],'r','linestyle','none','facealpha',0.3);
hp2 = patch([xo2,xo1],[yo2,yo1],'g','linestyle','none','facealpha',0.3);

hold off
saveas(figureOrbits,'lgs-orbit-comparison.png')