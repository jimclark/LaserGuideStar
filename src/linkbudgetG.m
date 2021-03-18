%% Gaussian beam link budget
% Assumption: Gaussian beam waist is 1/3 of the Tx aperture radius (Dtx/6)
% Provide all inputs in SI units (W, m, m, m, m).
% BWtx is FWHM beamwidth
% Includes 3 dB pointing loss (i.e. +/- BWtx/3.4)
function [Prx,Photrx,appMag,BWtx] = linkbudgetG(Ptx,Dtx,range,lambda,Drx)

c = 299792458;
h = 6.626e-34;

w0 = Dtx/6;
zR = pi.*(w0.^2)./lambda;

wZ = w0.*sqrt(1+(range./zR).^2);

BWtx = (2*atan(wZ./range))/1.7; % will converge to (2*lambda/(pi*w0))/1.7
fluxrec = 2*Ptx./(pi.*wZ.^2);

Arx = 0.25.*pi.*Drx.^2;
Prx = (10^-0.3).*Arx.*fluxrec; % 3 dB pointing loss accounted here.

Ephot = h.*c./lambda;
Photrx = Prx./Ephot;

Fx0 = zeros(size(lambda));
BWwideband = zeros(size(lambda));

for i = 1:numel(lambda)
    
    if (lambda(i) < 398e-9 && lambda(i) > 332e-9)     % U
        Fx0(i) = 1810;
        BWwideband(i) = 66/365;
    elseif (lambda(i) < 492e-9)                    % B
        Fx0(i) = 4260;
        BWwideband(i) = 94/445;
    elseif (lambda(i) < 595e-9)                     % V
        Fx0(i) = 3640;
        BWwideband(i) = 88/551;
    elseif (lambda(i) < 727e-9)                    % R
        Fx0(i) = 3080;
        BWwideband(i) = 138/658;
    elseif (lambda(i) < 880.5e-9)                  % I
        Fx0(i) = 2550;
        BWwideband(i) = 149/806;
    elseif (lambda(i) < 1080e-9)                   % Y, e.g. Starshot
        Fx0(i) = 2075; % ESTIMATE!!!!
        BWwideband(i) = 120/1020;
    elseif (lambda(i) < 1326.5e-9)                 % J, e.g. ABL
        Fx0(i) = 1600;
        BWwideband(i) = 213/1220;
    elseif (lambda(i) < 1783.5e-9)                 % H, e.g. NODE/FLARE
        Fx0(i) = 1080;
        BWwideband(i) = 307/1630;
    elseif (lambda(i) < 2385e-9)                   % K
        Fx0(i) = 670;
        BWwideband(i) = 390/2190;
    end
    
end

FJy = Photrx.*h.*1e26./(BWwideband.*Arx);

appMag = -2.5.*log10(FJy./Fx0);

end