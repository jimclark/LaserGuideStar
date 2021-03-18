%% Photon noise calcs

LGS_size = 0.3; % 30 cm height
LGS_ang_size = LGS_size/range_LGS; % 7 nrad = 1.5 mas
planet_ang_size = 2*Re/(10*parsec); % 40 picorad = 8 uas
worst_case_planet_ang_size = 2*11*Re/(4.3*ly); % Jupiter around Alpha Cen would be 3 nrad (0.6 mas) across, still smaller than LGS.

area = 5*0.01*0.01; % 1 cm2 of area
flux = 1368; % W/m2
div = deg2rad(0.5); % divergence from solar size, flat reflector
scopearea = 0.25*pi*scope_d^2;
spotsize = range_LGS*div;
spotarea = pi*0.25*spotsize^2;
spherearea = 4*pi*0.25*range_LGS^2;
fluxrec = flux*area/spotarea;

recangle = scope_d/range_LGS;
recSA = 2*pi*(1-cos(recangle/2));

Imax = flux/pi;
fluxrec2 = Imax*recSA*area/scopearea;

fluxrecSilly = flux*area/spherearea; % What if we just radiate isotropically

areaFull = 0.2*0.3; % What if we put a full 6U reflection down LUVOIR?
fluxrecFull = flux*areaFull/spotarea;


flux0jy = 3640; % zero magnitude in V band, per http://www.astro.umd.edu/~ssm/ASTR620/mags.html, janskys
flux0pps = 1.51e7*flux0jy*0.16;
Ephot = h*c/(500e-9);
flux0wm2 = flux0pps*Ephot;
mag = -2.5*log10(fluxrec/flux0wm2);
mag2 = -2.5*log10(fluxrec2/flux0wm2);
magSilly = -2.5*log10(fluxrecSilly/flux0wm2);

magFull = -2.5*log10(fluxrecFull/flux0wm2);

% IR emissions, sum up to 2500 nm (LUVOIR sensitivity)

lambdas = (100:2500)*1e-9;
Blams = PlancksLaw(lambdas,273+30);
pow_ir = 0.03*(sum(Blams)*1e-9)*4*pi*(4*.2*.3); % W/m^2/sr/nm * nm * sr * m^2 = W
% Summing over all lambda, but really it's just 1700 nm+ that matters (i.e.
% K band).  The factor of 0.03 is the IR emissivity of aluminum.
pow_ir_rec = pow_ir*(0.25*pi*scope_d^2)/(4*pi*range_LGS^2); % how much IR received

Fx0 = 670;
BWwideband = 390/2190;
Arx = 0.25*pi*scope_d^2;

Ephot = h.*c./2385e-9;
Photrx = pow_ir_rec./Ephot;

FJy = Photrx.*h.*1e26./(BWwideband.*Arx);

appMag = -2.5.*log10(FJy./Fx0);