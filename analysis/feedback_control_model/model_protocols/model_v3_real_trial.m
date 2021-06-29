function [swim,grspeed,brain_state] = model_v3_real_trial (par, dt, reaf, rest_dur)
% sensory processing delay = 0.22 s
trial_dur=15;
num_frames=round((trial_dur+2*rest_dur)/dt);
grspeed=zeros(1,num_frames);
gr_starts=round(rest_dur/dt);
gr_ends=round((rest_dur+trial_dur)/dt);
grspeed(gr_starts:gr_ends)=10;
swim=false(1,num_frames);
brain_state=zeros(5,num_frames);
bout_counter=0;
frames_since_bout_start=0;
frames_since_bout_end=inf;
reaf0 = [1 0 0 0 0];
nob = size(reaf,1);
reaf=[reaf0; reaf];
sens_delay=round(0.22/dt);
for t=sens_delay+1:num_frames
    if bout_counter<=nob
        this_reaf=reaf(bout_counter+1,:);
    else
        this_reaf=reaf0;
    end
    [grspeed(t), frames_since_bout_start, frames_since_bout_end] = exp_iteration_v3 (grspeed(t), swim(t-1), frames_since_bout_start, frames_since_bout_end, this_reaf);
    [swim(t), brain_state(:,t)] = brain_iteration_v3(swim(t-1), brain_state(:,t-1), grspeed(t-sens_delay), par);
    if swim(t)
        if ~swim(t-1)
            bout_counter=bout_counter+1;
        end
    end
end