close all


% Jim's first DRM
% num_stars = 350;
% total_obs = 450; % 350 + revisit top 100

% Jim's second DRM (actually the first computed using this code) based on hearsay (350 stars, 1000 observations)
% num_stars = 350;
% total_obs = num_pts*2.8; % 2 observations for all targets, plus 8 more (10 total) for top 10%


exp_time = mean(obs_dur)*daysec; % Average time from Chris's schedule
exp_acc = max_bg_acc; % From OrbitCalcs2dome3 calculations

max_simult_lgs = 60; % Let's say we want as many as 30 of these things active at once...

[lgs_req_std,domains_req_std,time_req_std] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass,sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);

[lgs_req_opt,domains_req_opt,time_req_opt] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass_opt_tot,sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);


sc_masses = 12:0.1:30;
lgs_req_mass = zeros(size(sc_masses));
domains_req_mass = zeros(size(sc_masses));
time_req_mass = zeros(size(sc_masses));

for i = 1:numel(sc_masses)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,range_LGS,sc_masses(i),sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req_mass(i) = lgs_req_temp;
    domains_req_mass(i) = domains_req_temp;
    time_req_mass(i) = time_req_temp;
end

sc_fuels = 1.0:0.1:3.0;

lgs_req_fuel = zeros(size(sc_fuels));
domains_req_fuel = zeros(size(sc_fuels));
time_req_fuel = zeros(size(sc_fuels));

for i = 1:numel(sc_fuels)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass_opt + sc_prop_dry_nom + sc_fuels(i),sc_fuels(i),sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req_fuel(i) = lgs_req_temp;
    domains_req_fuel(i) = domains_req_temp;
    time_req_fuel(i) = time_req_temp;
end

ranges_LGS = 10e6:1e6:100e6;

lgs_req_range = zeros(size(ranges_LGS));
domains_req_range = zeros(size(ranges_LGS));
time_req_range = zeros(size(ranges_LGS));

for i = 1:numel(ranges_LGS)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,ranges_LGS(i),sc_mass_opt_tot,sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req_range(i) = lgs_req_temp;
    domains_req_range(i) = domains_req_temp;
    time_req_range(i) = time_req_temp;
end

max_simult_lgses = 1:30;

lgs_req_simult = zeros(size(max_simult_lgses));
domains_req_simult = zeros(size(max_simult_lgses));
time_req_simult = zeros(size(max_simult_lgses));

for i = 1:numel(max_simult_lgses)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass_opt_tot,sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgses(i));
    lgs_req_simult(i) = lgs_req_temp;
    domains_req_simult(i) = domains_req_temp;
    time_req_simult(i) = time_req_temp;
end

% Looking at increasing the number of stars, observations, and mission
% times...getting some weird results here.  Going from 19 needed for stock
% LUVOIR to 43 needed for observing twice as many stars (and twice as many
% observations) in the same time is reasonable, but it's saying we can get
% that done in 1.9 years when we needed 3.5 for stock LUVOIR is a little
% weird...something for future work.

nums_stars = [num_stars, num_stars, num_stars*2, num_stars*2];
totals_obs = [total_obs, total_obs*2, total_obs*2, total_obs*2];
total_mission_times = [total_mission_time, total_mission_time*2, total_mission_time, total_mission_time*2];

new_missions_titles = categorical({'1x stars, 1x obs, 1x dur (stock)','1x stars, 2x obs, 2x dur','2x stars, 2x obs, 1x dur','2x stars, 2x obs, 2x dur'});

lgs_req_nstars = zeros(size(nums_stars));
domains_req_nstars = zeros(size(nums_stars));
time_req_nstars = zeros(size(nums_stars));

for i = 1:numel(nums_stars)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(nums_stars(i),totals_obs(i),range_LGS,sc_mass_opt_tot,sc_fuel_nom,sc_isp_nom,sc_max_thrust_nom,total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req_nstars(i) = lgs_req_temp;
    domains_req_nstars(i) = domains_req_temp;
    time_req_nstars(i) = time_req_temp;
end

figuremass = figure;
plot(sc_masses,lgs_req_mass,sc_mass_opt_tot,lgs_req_opt,'r*','linewidth',2)
legend('Total mass trade','Baseline case')
title('Minimum LGS required vs. LGS mass')
xlabel('LGS mass (kg)')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figuremass,sprintf('DRM-sensit-mass-%.2gk.png',range_LGS/1e6))

figurefuel = figure;
plot(sc_fuels,lgs_req_fuel,sc_fuel_nom,lgs_req_opt,'r*','linewidth',2)
legend('Fuel mass trade','Baseline case')
title('Minimum LGS required vs. LGS fuel mass')
xlabel('LGS fuel mass (kg)')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figurefuel,sprintf('DRM-sensit-fuel-%.2gk.png',range_LGS/1e6))

figurerange = figure;
plot(ranges_LGS/1e6,lgs_req_range,range_LGS/1e6,lgs_req_opt,'r*','linewidth',2)
legend('LGS range trade','Baseline case')
title('Minimum LGS required vs. Telescope-LGS range')
xlabel('Telescope-LGS range (1000''s km)')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figurerange,sprintf('DRM-sensit-range-%.2gk.png',range_LGS/1e6))

figuresimult = figure;
plot(max_simult_lgses,lgs_req_simult,domains_req_opt,lgs_req_opt,'r*','linewidth',2)
legend('Max active trade','Baseline case')
title('Minimum LGS required vs. Max LGS simult.')
xlabel('Maximum LGS simultaneously active')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figuresimult,sprintf('DRM-sensit-domains-%.2gk.png',range_LGS/1e6))

figuresimulttime = figure;
plot(max_simult_lgses,time_req_simult/yrsec,domains_req_opt,time_req_opt/yrsec,'r*',max_simult_lgses,5*ones(size(max_simult_lgses)),'linewidth',2)
legend('Max active trade','Baseline case','5-year limit')
ylim([0 17])
title('Time required for campaign vs. Max LGS simult.')
xlabel('Maximum LGS simultaneously active')
ylabel('Years to execute survey campaign')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figuresimulttime,sprintf('DRM-sensit-domains-time-%.2gk.png',range_LGS/1e6))

figurenstars = figure;
bar(new_missions_titles,lgs_req_nstars)