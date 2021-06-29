function [grspeed, frames_since_bout_start, frames_since_bout_end] = exp_iteration_v3 (grspeed, swim, frames_since_bout_start, frames_since_bout_end, reaf)
% dt = 0.005
% reaf = [gain, lag(frames), shunted(1/0), gain drop start, gain drop end]
if swim
    frames_since_bout_start=frames_since_bout_start+1;
    frames_since_bout_end=0;
    if frames_since_bout_start>reaf(2) % if bout is longer than lag
        if frames_since_bout_start<reaf(4) || frames_since_bout_start>reaf(5) % if gain is not "dropped"
            grspeed=grspeed-reaf(1)*20;
        end
    end
else
    frames_since_bout_end=frames_since_bout_end+1;
    frames_since_bout_start=0;
    if reaf(3)==0
        if frames_since_bout_end<=reaf(2)
            grspeed=grspeed-reaf(1)*20;
        end
    end
end
