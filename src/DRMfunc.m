function [lgs_req,domains_req,time_req] = DRMfunc(num_stars,total_obs,range,sc_mass,sc_fuel,sc_isp,sc_max_thrust,total_mission_time,exp_time,exp_acc,max_simult_lgs)

g0 = 9.8066;

idxs = [1 3]; % Only looking at the first and third "stars" in the Fibonacci spiral
phi = acos(1-2.*(idxs-0.5)./num_stars);
theta = pi*(1+sqrt(5))*(idxs-0.5);
x = cos(theta).*sin(phi);
y = sin(theta).*sin(phi);
z = cos(phi);

typSep = sqrt((x(1)-x(2)).^2+(y(1)-y(2)).^2+(z(1)-z(2)).^2);

sc_dv = sc_isp*g0*log(sc_mass/(sc_mass-sc_fuel));
sc_max_acc = sc_max_thrust/sc_mass;

min_maneuver_time = 2*sqrt(range*typSep/sc_max_acc);

desired_obs_interval = total_mission_time/total_obs;
exp_dv = exp_time*exp_acc;

total_req_lgs_cons = zeros(1,max_simult_lgs); % conservative
total_req_lgs_opt = zeros(1,max_simult_lgs); % optimistic
total_req_time = zeros(1,max_simult_lgs);
nums_domains = 1:max_simult_lgs;

for i = 1:max_simult_lgs
    num_domains = nums_domains(i);
    exp_per_domain = total_obs/num_domains;
    min_exp_interval = min_maneuver_time/num_domains;
    
    speed_factor = (desired_obs_interval-exp_time)/min_exp_interval;
    
    if speed_factor > 1
        maneuver_time = min_maneuver_time*speed_factor;
        sc_acc = sc_max_acc/(speed_factor^2);
    else
        maneuver_time = min_maneuver_time;
        sc_acc = sc_max_acc;
    end
    
    total_exp_dv = maneuver_time*sc_acc + exp_dv;
    
    exp_per_lgs = floor(sc_dv/total_exp_dv);
    
    total_req_lgs_cons(num_domains) = num_domains*ceil(exp_per_domain/exp_per_lgs);
    total_req_lgs_opt(num_domains) = max(ceil(num_domains*exp_per_domain/exp_per_lgs),num_domains);
    total_req_time(num_domains) = exp_per_domain*(exp_time + maneuver_time);
end

lgs_req = min(total_req_lgs_cons);
domains_req = find(total_req_lgs_cons==min(total_req_lgs_cons),1,'last');
time_req = total_req_time(domains_req);

if (time_req > total_mission_time)
    [time_req,domains_req] = min(total_req_time);
    lgs_req = total_req_lgs_opt(domains_req);
end

end