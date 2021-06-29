function [swim,grspeed,brain_state] = model_one_bout_trial_v3 (par, dt)
running=true;
sens_delay=round(0.22/dt);
t=sens_delay;
t_from_start_to_gr=0.7/dt;
t_from_bout_end_to_gr_end=0.4/dt;
reaf0 = [1 0 0 0 0];
brain_state=zeros(5,t);
swim=false(1,t);
frames_since_bout_start=0;
frames_since_bout_end=0;
fish_swam=false;
while running
    t=t+1;
    if t>t_from_start_to_gr && (frames_since_bout_end<t_from_bout_end_to_gr_end || ~fish_swam)
        grspeed(t)=10;
    else
        grspeed(t)=0;
    end
    [grspeed(t), frames_since_bout_start, frames_since_bout_end] = exp_iteration_v3 (grspeed(t), swim(t-1), frames_since_bout_start, frames_since_bout_end, reaf0); 
    [swim(t), brain_state(:,t)] = brain_iteration_v3(swim(t-1), brain_state(:,t-1), grspeed(t-sens_delay), par);
    if swim(t-1)
        fish_swam=true;
    end
    if frames_since_bout_end > t_from_bout_end_to_gr_end + t_from_start_to_gr && fish_swam
        running=false;
    end
    if t>=10/dt
        running=false;
    end
end