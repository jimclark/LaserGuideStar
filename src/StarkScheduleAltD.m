% This one is based on the TSP single-salesman circuit -- each LGS will
% pick up stars until it can't anymore, and then hand off to the next.
% Trying to fit a full circuit into 1/6th the time and dV budget.

% Needs to run ham_StarkSkymap first.

close all

lgs_count = 1;
lgs_dv = dvcap/6;
lgs_time = 0;

lgs_dvs = dvcap/6;
lgs_times = 0;

exp_time = mean(obs_dur)*daysec; % Average time from Chris's schedule
exp_dv = avg_acc*exp_time;

obs_asgn = zeros(size(starids));

speed_factor = 0.19; % reduce this to force the LGS to transit more slowly than max speed
% Empirically works down to 13 sats required (yay!) at ~0.19

sc_max_acc = (sc_max_thrust_nom/sc_mass_opt_tot)*speed_factor;

for i = 1:numel(starids)
    
    
    if (exp_dv > lgs_dv) || (lgs_time + exp_time > total_mission_time/6)
        lgs_dvs(lgs_count) = lgs_dv;
        lgs_times(lgs_count) = lgs_time;
        lgs_count = lgs_count + 1;
        lgs_dv = dvcap/6;
        lgs_time = 0;
    end
    
    obs_asgn(i) = lgs_count;
    
    lgs_dv = lgs_dv - exp_dv;
    lgs_time = lgs_time + exp_time;
    
    curr_starid = starids(hamStark(i));
    curr_starlat = starlats(hamStark(i));
    curr_starlon = starlons(hamStark(i));
    
    next_starid = starids(hamStark(i+1));
    next_starlat = starlats(hamStark(i+1));
    next_starlon = starlons(hamStark(i+1));
    
    next_dist = deg2rad(distance(curr_starlat,curr_starlon,next_starlat,next_starlon));
    
    transit_time = 2*sqrt(range_LGS*next_dist/sc_max_acc);
    transit_dv = transit_time*sc_max_acc;
    
    if (transit_dv < lgs_dv) && (lgs_time + transit_time < total_mission_time/6)
        lgs_dv = lgs_dv - transit_dv;
        lgs_time = lgs_time + transit_time;
    else
        lgs_dvs(lgs_count) = lgs_dv;
        lgs_times(lgs_count) = lgs_time;
        lgs_count = lgs_count + 1;
        lgs_dv = dvcap/6;
        lgs_time = 0;
    end
end

lgs_dvs(lgs_count) = lgs_dv;
lgs_times(lgs_count) = lgs_time;

num_sats = max(obs_asgn);

disp(num_sats)

C_obs = cell(num_sats,3);

obs_per_sat = zeros(size(1:num_sats));

for i = 1:num_sats
    obs_per_sat(i) = sum(obs_asgn==i);
    
    obs_made = starids(hamStark(find(obs_asgn == i)));
    
    obs_lats = starlats(hamStark(find(obs_asgn == i)));
    obs_lons = starlons(hamStark(find(obs_asgn == i)));
    C_obs{i,1} = obs_made;
    C_obs{i,2} = obs_lats;
    C_obs{i,3} = obs_lons;
end

figOPS = figure;
plot(6*obs_per_sat,'linewidth',2)
ylim([0 6*25])
title('Nr. obs. supported by each LGS satellite')
xlabel('LGS satellite number')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figOPS,'StarkSchedule_tsp_obs_per_sat.png')

figCloud = figure;
plot(obs_asgn,'x','linewidth',2)
title('Observations supported by each LGS satellite')
xlabel('Observation number')
ylabel('LGS satellite number')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figCloud,'StarkSchedule_tsp_obs_cloud.png')

figTime = figure;
hold on
plot(6*lgs_times/daysec,'linewidth',2)
plot([0 num_sats],[total_mission_time/(daysec) total_mission_time/(daysec)],'linewidth',2)
hold off
ylim([0 365])
title('Days of engagement by LGS spacecraft')
legend('LGS operation time','Max mission duration')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figTime,'StarkSchedule_tsp_time_remaining.png')

figDV = figure;
hold on
plot(lgs_dvs*6,'linewidth',2)
plot([0 num_sats],[dvcap dvcap],'linewidth',2)
hold off
title('dV remaining in each LGS spacecraft')
legend('LGS dV remaining','Initial dV capacity')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figDV,'StarkSchedule_tsp_dv_remaining.png')


colors_for_plot = lines(num_sats);

figureMap = figure;
axesm('MapProjection','robinson','Grid','on','GLineWidth',2)
% axesm('MapProjection','stereo','MapLatLimit',[-83 -90],'PLineLocation',1,'ParallelLabel','on','Grid','on','GLineWidth',2)
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2);
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
% legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
for i = 9
    plotm(C_obs{i,2},C_obs{i,3},'Color',colors_for_plot(i,:), 'linewidth', 2)
    
end

set(gca, 'fontsize', 14,'linewidth',2)
saveas(figureMap,'StarkSchedule_tsp_map.png')

figureGlobe = figure;

axesm('globe','Grid','on','GLineWidth',2,'MeridianLabel','on','MLabelParallel','equator','ParallelLabel','on','PLabelMeridian','prime')
p1 = scatterm(starlats,starlons,'*', 'linewidth', 2);
p2 = scatterm(deeplats,deeplons,'rv', 'linewidth', 2);
p3 = scatterm(brightlats,brightlons,'g+', 'linewidth', 2);
% legend([p1 p3 p2],{'Stark 2015 targets','Magnitude 2 stars','Hubble/Chandra deep fields'})
set(gca, 'fontsize', 14,'linewidth',2)

for i = 1:num_sats
    plot3m(C_obs{i,2},C_obs{i,3},0.02*ones(size(C_obs{i,1})),'Color',colors_for_plot(i,:), 'linewidth', 2)
    
end

base = zeros(180,360);
baseR = georefcells([-90 90],[0 360],size(base));
copperColor = [0.62 0.38 0.24];
geoshow(base,baseR,'FaceColor',copperColor)
camlight right
material([.8 .9 .4])

saveas(figureGlobe,'StarkSchedule_tsp_globe.png')
