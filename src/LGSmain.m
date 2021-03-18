%% LGS Main Function

%% Setting up variables

AU = 1.496e11;
Re = 6371000;
g0 = 9.8066;
daysec = 60*60*24;
yrsec = daysec*365.25;
solarConst = 1366;
c = 299792458;
h = 6.626e-34;
ly = c*yrsec;
parsec = AU/(deg2rad(1/3600));

% Updated mission based on Chris Stark's figures.
num_stars = 259;
total_obs = 1539;
total_mission_time = 5*yrsec; % 5 years

% Telescope parameters
scope_d = 9.2; % meters
scope_seg_d = 1.15;
obs_lam = 500e-9; % visible observation

% Laser parameters
lambda = 980e-9; % staying out of the visible band
pwr_laser = 5;
D_laser = 0.035; % 3.5 cm aperture

% LGS parameters

sc_mass = 24; % Maximum mass of 12U CubeSat
sc_mass_opt = 11.5; % Mass of selected components, without propulsion system

% Propulsion system options
%prop_names = {'Accion TILE 5000 x2','Apollo Constellation Engine','Busek BIT-3','Enpulsion IFM Nano (typ. setting) x2','Enpulsion IFM Nano (max I_{sp}) x2','Phase Four x2','VACCO Green Monopropellant','VACCO Cold Gas'};
prop_names = {'Accion 5000','Apollo CE','Busek BIT-3','Enp. Nano (nom.)','Enp. Nano (max I_{sp})','Phase Four','VACCO Grn','VACCO Cold Gas'};
sc_fuel = [0.64, 1.0, 1.5, 0.46, 0.46, 1, 2, 1.03];
sc_prop_dry = [2.2, 4.5, 1.4, 1.44, 1.44, 3, 3, 2.46];
sc_isp = [1500, 1500, 2300, 3000, 6000, 900, 170, 75]; % Can flex Enpulsion up to 6000 sec by reducing thrust.
sc_max_thrust = 1e-3*[3, 33, 1.24, 0.7, 0.5, 2, 0.4, 0.1];

%% Derived parameters

range_LGS = scope_d^2/(2*lambda); % quarter-wave curvature across the telescope mirror -- 43,000 km

iwa_box_rad = range_LGS*(0.25*obs_lam/scope_d); % Radius of "target box" trying to stay waaaay inside coronagraph, 0.25 lambda/D
iwa_box_rad_relax = range_LGS*(obs_lam/scope_d); % Radius of "target box" trying to stay kinda inside coronagraph, 1 lambda/D
seg_box_rad = range_LGS*(lambda/scope_seg_d); % Radius of "target box" trying to keep wavefronts flat (i.e. no more than one wavelength error) on each mirror segment

div_laser = 2*(lambda/(pi*(D_laser/6))); % 107 urad (980 nm), full-width Gaussian divergence
% Gaussian beam waist is 1/3rd the actual diameter of the main optic
% Todo: break out gaussian beam divergence into its own calculation

[Prx,Photrx,appMag,bw] = linkbudgetG(pwr_laser,D_laser,range_LGS,lambda,scope_d);
[Prx2,Photrx2,appMag2,bw2] = linkbudgetG(pwr_laser,D_laser,scope_d^2/(2*532e-9),532e-9,scope_d);

idxs = [1 3]; % Only looking at the first and third "stars" in the Fibonacci spiral
phi = acos(1-2.*(idxs-0.5)./num_stars);
theta = pi*(1+sqrt(5))*(idxs-0.5);
x = cos(theta).*sin(phi);
y = sin(theta).*sin(phi);
z = cos(phi);

sep_std = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2+(z(1)-z(2)).^2); % standard unit-sphere separation between adjacent stars (scale by range)

%% Preliminary propulsion-system selection

single_maneuver_time = 2.*sqrt(range_LGS*sep_std*sc_mass./sc_max_thrust);
smt_opt = 2.*sqrt(range_LGS*sep_std*(sc_mass_opt+sc_prop_dry+sc_fuel)./sc_max_thrust);
smt_opt_days = smt_opt/(60*60*24);
smt_opt_dv = 2.*sqrt(range_LGS*sep_std*sc_max_thrust./(sc_mass_opt+sc_prop_dry+sc_fuel));
single_maneuver_days = single_maneuver_time./(60*60*24);
single_maneuver_dv = 2.*sqrt(range_LGS*sep_std.*sc_max_thrust./sc_mass);
dv_caps = g0.*sc_isp.*log(sc_mass./(sc_mass-sc_fuel));
dv_caps_opt = g0.*sc_isp.*log((sc_mass_opt+sc_prop_dry+sc_fuel)./(sc_mass_opt+sc_prop_dry));
number_maneuvers = floor(dv_caps./single_maneuver_dv);
num_man_opt = dv_caps_opt./smt_opt_dv;
slowest_ep_smt_opt_days = max(smt_opt_days.*(sc_isp>1000));
num_man_opt_et = num_man_opt.*(slowest_ep_smt_opt_days./smt_opt_days).*(smt_opt_days<=slowest_ep_smt_opt_days); % Scale by the maneuver time of the lowest-thrust electric prop system (i.e. the longest maneuver time with a prop system isp > 1000 sec).

[max_mans_opt,idx_max_mans_opt] = max(num_man_opt_et);

fprintf("Selected propulsion system: %s\n",prop_names{idx_max_mans_opt});

sc_prop_dry_nom = sc_prop_dry(idx_max_mans_opt);
sc_fuel_nom = sc_fuel(idx_max_mans_opt);
sc_isp_nom = sc_isp(idx_max_mans_opt);
sc_max_thrust_nom = sc_max_thrust(idx_max_mans_opt);

sc_mass_opt_tot = sc_mass_opt+sc_prop_dry_nom+sc_fuel_nom; % total opt mass
dvcap = g0*sc_isp_nom*log(sc_mass_opt_tot/(sc_mass_opt_tot-sc_fuel_nom));

%% Run other scripts

OrbitCalcs2dome3 % to get the maximum background acceleration
OrbitCalcs2dome2
OrbitCalcs2dome
NoiseCalcsPropSens
NoiseCalcs
StarkSkymap

%%
DRM_prop_options
DRM_sensitivity

%%

if isfile('stark_skymap_tsp.mat')
    load('stark_skymap_tsp.mat')
else
    StarkSkymap_TSP % This can take some time.
end

if isfile('stark_skymap_tsp_ham.mat')
    load('stark_skymap_tsp_ham.mat')
else
    ham_StarkSkymap
end

%%

StarkSchedule
StarkScheduleAltB
StarkScheduleAltD

%%
SkyCalcs
SkyCalcsOffGEO
SkyCalcsOffGEO2

%%
HEO_LGS
SkyCalcsHEO
PowerCalcs
OrbitCalcs3CL3
LGSretro