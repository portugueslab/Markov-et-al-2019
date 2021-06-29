function [bd, id] = model_compute_parameters_v3 (swim,dt)
bout_starts_ids=find(diff([0 swim])==1);
bout_ends_ids=find(diff([0 swim 0])==-1);
if length(bout_starts_ids)==3 && length(bout_ends_ids)==3
    bd=(bout_ends_ids(2)-bout_starts_ids(2))*dt;
    id=(bout_starts_ids(3)-bout_ends_ids(2))*dt;
else
    bd=nan;
    id=nan;
end