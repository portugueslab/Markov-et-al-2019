%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_PC_imaging_results=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_trigaver=[pathname_PC_imaging_results 'triggered_traces\'];
pathname_behavior=[pathname_PC_imaging_results 'behavior\'];
pathname_scores=[pathname_PC_imaging_results 'scores\'];

%% loop through all fish
cd(pathname_trigaver);
all_fish=dm_dir('*_f*_trig_traces.mat');
all_fish=strrep(all_fish,'_trig_traces.mat','');
n_fish=length(all_fish);
fish_id = all_fish{1};
load([pathname_trigaver fish_id '_trig_traces.mat']);
dt=time_trig(2)-time_trig(1);
dt=round(dt*100)/100;
clipping=1.2;
progressbar('Computing scores...');
for f=1:n_fish
    fish_id = all_fish{f};
    filename = [pathname_scores fish_id '_scores.mat'];
    load([pathname_trigaver fish_id '_trig_traces.mat']);
    scores_im = compute_scores_here(traces_bout_trig,time_trig,clipping);
    scores_motor_regr = compute_scores_here(trace_motor_regr_bout_trig,time_trig,clipping)';
    load([pathname_behavior fish_id '_behavior.mat'],'trials');
    scores_be=trials.bout_duration.first; 
    save(filename,'scores_be','scores_im','scores_motor_regr')
    progressbar(f/n_fish);
end

function [scores] = compute_scores_here(data,time,clipping)
scores=squeeze(nanmean(data(:,time>0 & time<=clipping,:),2));
end