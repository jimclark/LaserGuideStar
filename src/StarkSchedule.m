% Ingests the list of observations from Stark, verbatim, rigidly fits them
% into a schedule of evenly-spaced observations over four years, and
% schedules LGSs accordingly.  Levels out around 34 LGSs required, should
% only need 20 or so.

% start by running StarkSkymap

close all;

obs_asgn = zeros(size(obs_cts));

lgs_count = 1;
lgs_dv = dvcap;

lgs_dvs = dvcap;

speed_factor = 0.23; % reduce this to force the LGS to transit more slowly than max speed
% Empirically, 0.23 seems to have the best results.

sc_max_acc = (sc_max_thrust_nom/sc_mass_opt_tot)*speed_factor;

while sum(obs_asgn==0) > 0
    
    unasgn = find(obs_asgn==0);
    
%     disp(numel(unasgn))
%     disp(lgs_count)
    
    curr_idx = unasgn(1);
    
    exp_dv = max_bg_acc*daysec*obs_dur(curr_idx);
    
    if(exp_dv > lgs_dv)
        lgs_dvs(lgs_count) = lgs_dv;
        lgs_count = lgs_count + 1;
        lgs_dv = dvcap;
    end
    
    lgs_dv = lgs_dv - exp_dv;
    obs_asgn(curr_idx) = lgs_count;
    
    curr_starid = obs_ids(curr_idx);
    curr_starlat = starlats(find(starids==curr_starid));
    curr_starlon = starlons(find(starids==curr_starid));
    
    next_obs = 2;
    
    while next_obs <= numel(unasgn)
        test_idx = unasgn(next_obs);
        test_starid = obs_ids(test_idx);
        
        test_starlat = starlats(find(starids==test_starid));
        test_starlon = starlons(find(starids==test_starid));
        
        test_dist = deg2rad(distance(curr_starlat,curr_starlon,test_starlat,test_starlon));
        
        transit_time_req = ((test_idx-curr_idx)*total_mission_time/total_obs)-daysec*obs_dur(curr_idx);
        transit_time_actual = 2*sqrt(range_LGS*test_dist/sc_max_acc);
        transit_dv = transit_time_actual*sc_max_acc + max_bg_acc*daysec*obs_dur(curr_idx);
        
        if (transit_time_actual < transit_time_req) && (transit_dv < lgs_dv)
            curr_idx = test_idx;
            lgs_dv = lgs_dv - transit_dv;
            obs_asgn(curr_idx) = lgs_count;
            
            curr_starid = obs_ids(curr_idx);
            curr_starlat = starlats(find(starids==curr_starid));
            curr_starlon = starlons(find(starids==curr_starid));
        end
        next_obs = next_obs + 1;
    end
    
    
    if sum(obs_asgn==0) > 0
        lgs_dvs(lgs_count) = lgs_dv;
        lgs_count = lgs_count + 1;
        lgs_dv = dvcap;
    end
    
end

lgs_dvs(lgs_count) = lgs_dv;

num_sats = max(obs_asgn);

disp(num_sats)

C_obs = cell(num_sats,3);

obs_per_sat = zeros(size(1:num_sats));

for i = 1:num_sats
    obs_per_sat(i) = sum(obs_asgn==i);
    
    obs_made = obs_ids(find(obs_asgn == i));
    obs_idxs = zeros(size(obs_made));
    for j = 1:obs_per_sat(i)
        obs_idxs(j) = find(starids == obs_made(j));
    end
    obs_lats = starlats(obs_idxs);
    obs_lons = starlons(obs_idxs);
    C_obs{i,1} = obs_idxs;
    C_obs{i,2} = obs_lats;
    C_obs{i,3} = obs_lons;
end

figOPS = figure;
plot(obs_per_sat,'linewidth',2)
title('Nr. obs. supported by each LGS satellite')
xlabel('LGS satellite number')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figOPS,'StarkSchedule_rigid_obs_per_sat.png')

figCloud = figure;
plot(obs_asgn,'x','linewidth',2)
title('Observations supported by each LGS satellite')
xlabel('Observation number')
ylabel('LGS satellite number')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figCloud,'StarkSchedule_rigid_obs_cloud.png')

% figure;
% hold on
% plot(lgs_times/daysec)
% plot([0 num_sats],[total_mission_time/daysec total_mission_time/daysec])
% hold off
% title('Days of engagement by LGS spacecraft')
% legend('LGS operation time','Max mission duration')

figDV = figure;
hold on
plot(lgs_dvs,'linewidth',2)
plot([0 num_sats],[dvcap dvcap],'linewidth',2)
hold off
title('dV remaining in each LGS spacecraft')
legend('LGS dV remaining','Initial dV capacity')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figDV,'StarkSchedule_rigid_dv_remaining.png')


colors_for_plot = lines(num_sats);

figureMap = figure;
axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
% axesm('MapProjection','stereo','MapLatLimit',[-83 -90],'PLineLocation',1,'ParallelLabel','on','Grid','on','GLineWidth',2)
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2);
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
% legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
for i = [1 10 floor(num_sats/7)*7]
    plotm(C_obs{i,2},C_obs{i,3},'Color',colors_for_plot(i,:), 'linewidth', 2)
    
end

set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,'StarkSchedule_rigid_map.png')
