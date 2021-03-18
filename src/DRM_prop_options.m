close all

% Eventual todo: use this script to pick which propulsion system goes
% into the "nom" values.

exp_time = mean(obs_dur)*daysec; % Average time from Chris's schedule
exp_acc = max_bg_acc; % From OrbitCalcs2dome3 calculations

max_simult_lgs = 100; % Let's say we want as many as 30 of these things active at once...

lgs_req = zeros(size(sc_fuel));
domains_req = zeros(size(sc_fuel));
time_req = zeros(size(sc_fuel));

for i = 1:numel(sc_fuel)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass,sc_fuel(i),sc_isp(i),sc_max_thrust(i),total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req(i) = lgs_req_temp;
    domains_req(i) = domains_req_temp;
    time_req(i) = time_req_temp;
end


lgs_req_opt = zeros(size(sc_fuel));
domains_req_opt = zeros(size(sc_fuel));
time_req_opt = zeros(size(sc_fuel));

for i = 1:numel(sc_fuel)
    [lgs_req_temp,domains_req_temp,time_req_temp] = DRMfunc(num_stars,total_obs,range_LGS,sc_mass_opt+sc_prop_dry(i)+sc_fuel(i),sc_fuel(i),sc_isp(i),sc_max_thrust(i),total_mission_time,exp_time,exp_acc,max_simult_lgs);
    lgs_req_opt(i) = lgs_req_temp;
    domains_req_opt(i) = domains_req_temp;
    time_req_opt(i) = time_req_temp;
end

prop_labels = categorical(prop_names);

figBar = figure;
bar(prop_labels,[lgs_req;lgs_req_opt]','linewidth',2)
legend('24 kg','11.5 kg + prop','Location','northwest')
title(sprintf('Propulsion trade, %.2g,000 km range',range_LGS/1e6))
ylabel('Min. LGSs required')
set(gca, 'fontsize', 14,'linewidth',2)

%%
saveas(figBar,sprintf('DRM_prop_options_%.2gk.png',range_LGS/1e6))
ylim([0 60])
saveas(figBar,sprintf('DRM_prop_options_%.2gk_zoom.png',range_LGS/1e6))