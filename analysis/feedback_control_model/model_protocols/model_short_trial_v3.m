function [swim,grspeed,brain_state] = model_short_trial_v3 (par, dt, reaf)
% sensory processing delay = 0.22 s
running=true;
sens_delay=round(0.22/dt);
gr_starts_frame=round(0.3/dt);
t=sens_delay;
reaf0 = [1 0 0 0 0; 1 0 0 0 0; reaf; 1 0 0 0 0];
num_frames=round(10/dt);
brain_state=zeros(5,num_frames);
swim=false(1,num_frames);
frames_since_bout_start=0;
frames_since_bout_end=inf;
bout_counter=0;
grspeed=zeros(1,num_frames);
grspeed(gr_starts_frame:end)=10;
while running
    t=t+1;
    [grspeed(t), frames_since_bout_start, frames_since_bout_end] = exp_iteration_v3 (grspeed(t), swim(t-1), frames_since_bout_start, frames_since_bout_end, reaf0(bout_counter+1,:));
    [swim(t), brain_state(:,t)] = brain_iteration_v3(swim(t-1), brain_state(:,t-1), grspeed(t-sens_delay), par);
    if swim(t)
        if ~swim(t-1)
            bout_counter=bout_counter+1;
            if bout_counter==3
                running=false;
            end
        end
    end
    if t>=num_frames
        running=false;
    end
end