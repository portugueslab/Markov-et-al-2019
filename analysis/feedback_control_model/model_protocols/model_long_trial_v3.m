function [swim,grspeed,brain_state] = model_long_trial_v3(par, dt, reaf)
% sensory processing delay = 0.22 s
num_frames=round(20/dt); % 20 sec
gr_starts=round(2.5/dt); % 2.5 sec
gr_ends=round(17.5/dt); % 17.5 sec
grspeed=zeros(1,num_frames);
grspeed(gr_starts:gr_ends)=10;
brain_state=zeros(5,num_frames);
swim=false(1,num_frames);
frames_since_bout_start=0;
frames_since_bout_end=inf;
sens_delay=round(0.22/dt);
for t=sens_delay+1:num_frames
    [grspeed(t), frames_since_bout_start, frames_since_bout_end] = exp_iteration_v3 (grspeed(t), swim(t-1), frames_since_bout_start, frames_since_bout_end, reaf);
    [swim(t), brain_state(:,t)] = brain_iteration_v3(swim(t-1), brain_state(:,t-1), grspeed(t-sens_delay), par);
end