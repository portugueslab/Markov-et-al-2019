function [trig_traces, time_trig] = dm_compute_triggered_traces(...
    trig_list,...
    s_pre,...
    s_post,...
    dt,...
    time_im,...
    traces)

id_pre = s_pre/dt;
id_post = s_post/dt;
trig_length=(s_pre+s_post)/dt;
time_trig = -s_pre+dt:dt:s_post;
time_trig=round(time_trig*100)/100;
num_trig = length(trig_list);
num_traces=size(traces,1);
trig_traces = nan(num_traces,trig_length,num_trig,'single');
for i=1:num_trig
    this_trig=trig_list(i);
    if ~isnan(this_trig)
        [~, this_trig_id] = min(abs(time_im-this_trig));
        trig_traces(:,:,i) = traces(:,this_trig_id-id_pre+1:this_trig_id+id_post);
    end
end
trig_traces = trig_traces - nanmean(trig_traces(:,1:id_pre,:),2);
end