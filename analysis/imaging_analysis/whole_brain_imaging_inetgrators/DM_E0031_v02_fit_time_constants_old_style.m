%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_behavior=[pathname_imaging_results 'behavior\'];
pathname_traces=[pathname_imaging_results 'traces\'];
pathname_trigaver=[pathname_imaging_results 'triggered_traces\'];
pathname_clustering=[pathname_imaging_results 'clustering\'];
pathname_time_constants=[pathname_imaging_results 'time_constants\'];

%% all fish to be analyzed
all_fish=dm_dir([pathname_traces '*_f*_traces.mat']);
all_fish=strrep(all_fish,'_traces.mat','');
n_fish=length(all_fish);

%% build regressors for different time constants
% time of forward moving grating
load([pathname_trigaver all_fish{1} '_trig_traces.mat'],'time_trig');
n_frames = length(time_trig);
% model activity of a neuron which is active during trials of forward
% moving grating and silent otherwise
model_act = zeros(1,n_frames);
model_act(time_trig>0)=1;
% build calcium kernel
load([pathname_behavior all_fish{1} '_behavior.mat'],'meta');
dt = time_trig(2)-time_trig(1);
s_half_decay = genotype2tau(meta.genotype);
my_kernel=2.^(-(dt:dt:10)/s_half_decay);
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
    load([pathname_trigaver fish_id '_trig_traces.mat'],'traces_gr_trig_mean');
    n_ROIs=size(traces_gr_trig_mean,1);
    load([pathname_clustering fish_id '_clustering.mat'],'sensmot_clust');
    ids_sens = find(sensmot_clust==1);
    traces_sens = traces_gr_trig_mean(ids_sens,:);
    n_ROIs_sens = length(ids_sens);
           
    %% compute correlation values between real and modelled traces, find the best tau for each sensory ROI
    X = traces_sens - nanmean(traces_sens,2);
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