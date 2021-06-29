%% pathnames
clc; close all; clear all;
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_PC_imaging_results=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_behavior=[pathname_PC_imaging_results 'behavior\'];
pathname_traces=[pathname_PC_imaging_results 'traces\'];
pathname_trigaver=[pathname_PC_imaging_results 'triggered_traces\'];

%% timing stuff
s_pre = 1;
s_post = 9;
dt=0.2; % 200 ms, corresponds to 5 Hz.

%% loop through all fish
all_fish=dm_dir([pathname_traces '*_f*_traces.mat']);
all_fish=strrep(all_fish,'_traces.mat','');
n_fish=length(all_fish);
for f = 1:n_fish
    fish_id = all_fish{f};
    filename = [pathname_trigaver fish_id '_trig_traces.mat'];
    
    %% create empty structure
    trigaver=struct;
    
    %% load the data
    load([pathname_traces fish_id '_traces.mat'],'traces','time_offsets','trace_motor_regr');
    load([pathname_behavior fish_id '_behavior.mat'],'bouts','grmov','meta','time_be');
    [n_ROIs, n_frames] = size(traces);
    time_im_ori=1/meta.F:1/meta.F:n_frames/meta.F;
    time_im=0:dt:time_im_ori(end)+30-dt;
    
    %% find trigger times
    gr_starts=time_be(diff(grmov)==1);
    bout1_starts = nan(1,length(gr_starts));
    for t=1:length(gr_starts)
        bouts_in_this_trial = bouts.start(bouts.trial==t);
        if ~isempty(bouts_in_this_trial)
            bout1_starts(t)= bouts_in_this_trial(1);
        end
    end
    
    %% interpolate traces
    traces_interp = zeros(n_ROIs,length(time_im),'single');
    for t=1:n_ROIs
        traces_interp(t,:) = interp1(time_im_ori+time_offsets(t), traces(t,:), time_im);
    end
    trace_motor_regr_interp = zeros(1,length(time_im),'single');
    trace_motor_regr_interp(1,:) = interp1(time_im_ori, trace_motor_regr, time_im);
    
    %% grating_onset
%     [trigaver.gr_on.data, trigaver.gr_on.time] = dm_compute_triggered_traces(gr_starts,s_pre,s_post,dt,time_im,traces_interp);
%     [trigaver_motor_regr.gr_on.data, trigaver_motor_regr.gr_on.time] = dm_compute_triggered_traces(gr_starts,s_pre,s_post,dt,time_im,trace_motor_regr_interp);
    
    %% first bout onset
    [traces_bout_trig, time_trig] = dm_compute_triggered_traces(bout1_starts,s_pre,s_post,dt,time_im,traces_interp);
    trace_motor_regr_bout_trig = dm_compute_triggered_traces(bout1_starts,s_pre,s_post,dt,time_im,trace_motor_regr_interp);
    
    %% save
    save(filename,'traces_bout_trig','trace_motor_regr_bout_trig','time_trig');
end