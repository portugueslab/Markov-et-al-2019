%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_behavior=[pathname_imaging_results 'behavior\'];
pathname_traces=[pathname_imaging_results 'traces\'];
pathname_clustering=[pathname_imaging_results 'clustering\'];
pathname_time_constants=[pathname_imaging_results 'time_constants\'];

%% all fish to be analyzed
all_fish=dm_dir([pathname_traces '*_f*_traces.mat']);
all_fish=strrep(all_fish,'_traces.mat','');
n_fish=length(all_fish);

%% build regressors for different time constants
% time of forward moving grating
load([pathname_behavior all_fish{1} '_behavior.mat'],'meta','time_be');
forw_gr_starts = 120+7.5:30:time_be(end);
forw_gr_ends=forw_gr_starts+15;
dt=1/meta.fs;
time_im = dt:dt:time_be(end);
n_frames = length(time_im);
% convert times into indices
forw_gr_starts=round(forw_gr_starts/dt);
forw_gr_ends=round(forw_gr_ends/dt);
% model activity of a neuron which is active during trials of forward
% moving grating and silent otherwise
model_act = zeros(1,n_frames);
for i=1:length(forw_gr_starts)
    model_act(forw_gr_starts(i):forw_gr_ends(i))=1;
end
% build calcium kernel
s_half_decay = genotype2tau(meta.genotype);
my_kernel=2.^(-(time_im)/s_half_decay);
% time constants to probe
max_tau=10;
probing_taus=dt:dt:max_tau;
n_taus=length(probing_taus);
% compute leaky integratrion for a list of taus and
% convolve it with the calcium kernel
integrated_traces = zeros(n_taus,n_frames);
% if tau == dt, no integration
for i = 1:n_taus
    tau=probing_taus(i);
    tau=dt/tau;
    this_integrated_trace = zeros(1,n_frames);
    if tau<inf
        for t=2:n_frames
            this_integrated_trace(t)=tau*model_act(t)-(tau-1)*this_integrated_trace(t-1);
        end
    else
        this_integrated_trace = model_act;
    end
    this_integrated_trace=conv(this_integrated_trace,my_kernel);
    integrated_traces(i,:) = this_integrated_trace(1:n_frames);
end
integrated_traces_mean = mean(integrated_traces,2);
Y = integrated_traces - integrated_traces_mean;
Y_sqr = sum(Y.^2,2);
   
%% loop through all fish
progressbar('Fish progress...','Interpolating traces...','Computing corrcoeffs...');
for f=1:n_fish
    fish_id=all_fish{f};
    filename = [pathname_time_constants fish_id '_time_constants.mat'];
    
    %% load the data
    load([pathname_traces fish_id '_traces.mat'],'traces','time_offsets');
    n_ROIs=size(traces,1);
    load([pathname_clustering fish_id '_clustering.mat'],'sensmot_clust');
    load([pathname_behavior fish_id '_behavior.mat'],'meta');
    im_dt=1/meta.fs;
    time_im_ori=im_dt:im_dt:size(traces,2)*im_dt;
    ids_sens = find(sensmot_clust==1);
    traces_sens = traces(ids_sens,:);
    n_ROIs_sens = length(ids_sens);
       
    %% interpolate all traces to time_im
    traces_interp = nan(n_ROIs_sens,n_frames,'single');
    for i=1:n_ROIs_sens
        traces_interp(i,:)=interp1(time_im_ori+time_offsets(i),traces_sens(i,:),time_im);
        progressbar([],i/n_ROIs_sens,[]);
    end
    
    %% compute correlation values between real and modelled traces, find the best tau for each sensory ROI
    traces_interp_mean = nanmean(traces_interp,2);
    X = traces_interp - traces_interp_mean;
    X_sqr = sum(X.^2,2,'omitnan');
    corrcoeffs  = nan(n_ROIs_sens,n_taus);
    for i = 1:n_taus
        corrcoeffs(:,i) = sum(X.*Y(i,:),2,'omitnan')./...
            sqrt(X_sqr*Y_sqr(i));
        progressbar([],[],i/n_taus);
    end
    [~,best_corrcoeff_ids]= max(corrcoeffs,[],2);
    time_constants_sens = probing_taus(best_corrcoeff_ids)';
    time_constants = nan(n_ROIs,1);
    time_constants(ids_sens) = time_constants_sens;
    
    %% save
    save(filename,'time_constants');
    progressbar(f/n_fish,[],[]);
end