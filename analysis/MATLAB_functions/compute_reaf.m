function [grspeed] = compute_reaf (swim, vigor, reaf)
% dt = 0.005
% reaf = [gain, lag(frames), shunted(1/0), gain drop start, gain drop end]
n=length(swim);
grspeed=10*ones(1,n);
start_frame=reaf(2)+1;
bout_frames=0;
for t=start_frame:n
    if swim(t-reaf(2))
        grspeed(t)=10-28.5*vigor(t-reaf(2))*reaf(1);
    end
    if swim(t)
        bout_frames=bout_frames+1;
    else
        if reaf(3)==1
            grspeed(t)=10;
        end
    end
    if reaf(4)~=0 || reaf(5)~=0
        if bout_frames>=reaf(4) && bout_frames<=reaf(5)
            grspeed(t)=10;
        end
    end
    if t==n
        a=1;
    end
end