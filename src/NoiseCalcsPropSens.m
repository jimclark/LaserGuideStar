%% Sensor and thruster noise calcs

close all;

max_upd_int = 6000;

upd_int_vec = 10:10:max_upd_int; % update interval

max_zerocost_pos_err = 0.25*max_bg_acc*upd_int_vec.^2; % run OrbitCalcs2dome3 first
typ_zerocost_pos_err = 0.25*avg_acc*upd_int_vec.^2;

figMaxPosErr = figure;
semilogy(upd_int_vec/60,max_zerocost_pos_err,'-',upd_int_vec/60,typ_zerocost_pos_err,'-',[0 max_upd_int/60],[seg_box_rad seg_box_rad],'--',[0 max_upd_int/60],[iwa_box_rad_relax iwa_box_rad_relax],'--',[0 max_upd_int/60],[iwa_box_rad iwa_box_rad],'--','linewidth', 2)
title('Maximum location error without penalty')
ylabel('LGS velocity-axis error (m)')
xlabel('Update interval (minutes)')
legend('Worst-case observation','Average observation',sprintf('Goal: +/- %.2g m (flat waves on segments)',seg_box_rad),sprintf('Goal: +/- %.2g m (stay in IWA, 1 λ/D)',iwa_box_rad_relax),sprintf('Goal: +/- %.2g m (stay in deep IWA, 0.25 λ/D)',iwa_box_rad),'Location','southeast')
set(gca, 'fontsize', 14,'linewidth',2)
saveas(figMaxPosErr,'max_pos_err_nocost.png')

seg_upd_time = 2*sqrt(seg_box_rad/max_bg_acc);
iwa_upd_time = 2*sqrt(iwa_box_rad/max_bg_acc);
iwa_upd_time_relax = 2*sqrt(iwa_box_rad_relax/max_bg_acc);

seg_upd_time_typ = 2*sqrt(seg_box_rad/avg_acc);
iwa_upd_time_typ = 2*sqrt(iwa_box_rad/avg_acc);
iwa_upd_time_typ_relax = 2*sqrt(iwa_box_rad_relax/avg_acc);

max_thr_err = 0.25*max_bg_thrust;
max_ang_err = asind(0.25);