% this progam performs pre-processing of the traces and ROI
% of the long-term adaptation experiment during PC lightsheet imaging.

%% pathnames
clc; close all; clear all;
% path to raw data
pathname_raw_data='J:\_Shared\experiments\E0030_long_term_adaptation\v09_lightsheet\data_imaging_suite2p\';
pathname_data_traces=[pathname_raw_data 'traces_filtered\'];
pathname_data_ori_planes=[pathname_raw_data 'plane_idxs\'];
pathname_data_ROIs=[pathname_raw_data 'rois_filtered\'];
pathname_data_raw_anat=[pathname_raw_data 'anatomies\'];
% path to pre-processed data uploaded with the manuscript
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_PC_imaging_results=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_behavior=[pathname_PC_imaging_results 'behavior\'];
pathname_traces=[pathname_PC_imaging_results 'traces\'];
pathname_ROIs=[pathname_PC_imaging_results 'ROIs\'];

%% cutoff frequency for high-pass filtering
fc_baseline=1/300; % Hz, 5 min period

%% loop through all fish
all_fish=dm_dir([pathname_data_traces '*_f*_traces.h5']);
all_fish=strrep(all_fish,'_traces.h5','');
n_fish=length(all_fish);
progressbar('Fish...','Traces...')
for f=1:n_fish
    fish_id=all_fish{f};
    
    %% load traces, ROIs, ori planes and anatomies, as well as  metadata and swim for motor regressor
    load([pathname_behavior fish_id '_behavior.mat'],'meta','swim','time_be');
    traces = h5read([pathname_data_traces fish_id '_traces.h5'],'/traces')';
    traces=double(traces);
    ROIs = tiff2mat(pathname_data_ROIs, [fish_id '_rois.tif']);
    ori_planes = h5read([pathname_data_ori_planes fish_id '_plane_idxs.h5'],'/plane_idxs')';
    anat = tiff2mat(pathname_data_raw_anat, [fish_id '_anatomy.tif']);
    [n_ROIs, n_frames] = size(traces);
    sz = size(anat);
    
    %% high-pass filter for removing slow drift in the traces
    [bb_baseline,aa_baseline]=butter(3,2*fc_baseline/meta.F,'high');
    traces(traces==0)=nan;
    
    %% process this fish
    ROI_coord_ori_lin = cell(n_ROIs,1);
    ROI_coord_ori_vx_rcp = cell(n_ROIs,1);
    time_offsets=nan(n_ROIs,1);
    traces2=nan(n_ROIs,n_frames,'single');
    for i=1:n_ROIs
        ROI_coord_ori_lin{i} = uint32(find(ROIs==i));
        [rows, columns, planes] = ind2sub(sz,ROI_coord_ori_lin{i});
        ROI_coord_ori_vx_rcp{i} = uint16([rows columns planes]);
        plane=unique(planes);
        time_offsets(i)=meta.pulse_times(ori_planes(plane));
        traces2(i,:) = process_trace_here(traces(i,:),bb_baseline,aa_baseline);
        progressbar([],i/n_ROIs)
    end
    traces=traces2; 
    
    %% build motor regressor
    time_im=1/meta.F:1/meta.F:n_frames/meta.F;
    s_half_decay = genotype2tau('GCaMP6s');
    my_kernel=2.^(-time_be/s_half_decay);
    trace_motor_regr=trace2regressor(swim,my_kernel,time_be,time_im);
    trace_motor_regr = process_trace_here(trace_motor_regr,bb_baseline,aa_baseline);
    
    %% save
    save([pathname_ROIs fish_id '_ROIs.mat'],'ROIs','ROI_coord_ori_lin','ROI_coord_ori_vx_rcp');
    save([pathname_traces fish_id '_traces.mat'],'traces','time_offsets','trace_motor_regr');
    progressbar(f/n_fish,[]);
end

function [trace] = process_trace_here(trace,bb_baseline,aa_baseline)
my_nans=isnan(trace);
trace(my_nans)=nanmean(trace);
trace=filtfilt(bb_baseline,aa_baseline,trace);
trace(my_nans)=nan;
trace=trace-nanmean(trace);
trace=trace/nanstd(trace);
trace=single(trace);
end