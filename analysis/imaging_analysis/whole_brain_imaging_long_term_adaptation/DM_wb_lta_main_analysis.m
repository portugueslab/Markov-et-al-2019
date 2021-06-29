%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data_code\data\';
pathname_wb_imaging_results=[pathname_MS_data 'whole_brain_imaging_long_term_adaptation\'];
pathname_behavior=[pathname_wb_imaging_results 'behavior\'];
pathname_traces=[pathname_wb_imaging_results 'traces\'];
pathname_processed_data=[pathname_wb_imaging_results 'processed_data\'];

%% timing stuff
s_pre = 1;
s_post = 4;
dt=0.2; % 200 ms, corresponds to 5 Hz.
time_trig=round((-s_pre+dt:dt:s_post)*100)/100;
n_frames_trig = length(time_trig);
score_tf=[time_trig>0;time_trig>0 & time_trig<2]; % sens tf; mot tf
baseline_tf=1:s_pre/dt-1;

%% parameters for bootstrapping
n_boots = 1000;
chunk_length = 23; % seconds
chunk_length = round(chunk_length/dt); % frames

%% for computing scores used for criteria
clipping=1.2;

%% loop through all fish
all_fish=dm_dir([pathname_traces '*_f*_traces.mat']);
all_fish=strrep(all_fish,'_traces.mat','');
n_fish=length(all_fish);
progressbar('Fish progress...','Interpolating traces...','Bootstrapping...','Fitting time constants...');
for f = 1:n_fish
    fish_id = all_fish{f};
    filename = [pathname_processed_data fish_id '_processed_data.mat'];
   
    %% load the behavioral data
    load([pathname_behavior fish_id '_behavior.mat'],'bouts','grmov','meta','time_be','trials');   
    
    %% load and interpolate traces
    load([pathname_traces fish_id '_traces.mat'],'traces','time_offsets','trace_motor_regr');
    [n_ROIs, n_frames_ori] = size(traces);
    im_dt = 1/meta.F;
    time_im_ori=(im_dt:im_dt:n_frames_ori/meta.F)-im_dt;
    time_im=dt:dt:time_im_ori(end)+30-dt;
    n_frames=length(time_im);
    traces_interp = zeros(n_ROIs,length(time_im),'single');
    for i=1:n_ROIs
        traces_interp(i,:) = interp1(time_im_ori+time_offsets(i), traces(i,:), time_im);
        progressbar([],i/n_ROIs,[],[])
    end
    trace_motor_regr_interp = zeros(1,length(time_im),'single');
    trace_motor_regr_interp(1,:) = interp1(time_im_ori, trace_motor_regr, time_im);
    clear traces trace_motor_regr;
      
    %% build triggers
    gr_starts=time_be(diff(grmov)==1);
    gr_ends=gr_starts+15;
    bout_starts=bouts.start;
    bout1_starts = nan(1,length(gr_starts));
    for t=1:length(gr_starts)
        bouts_in_this_trial = bout_starts(bouts.trial==t);
        if ~isempty(bouts_in_this_trial)
            bout1_starts(t)= bouts_in_this_trial(1);
        end
    end
    % build all triggers together
    all_trig=sort([gr_starts, gr_ends, bout_starts]);
    all_trig_cell = {gr_starts, bout_starts};
    num_trig_types=length(all_trig_cell);
    % work with triggers
    n_trig=zeros(1,num_trig_types);
    id_start=cell(1,num_trig_types);
    id_end=id_start;
    trig_length=id_start;
    for i=1:num_trig_types
        % remove triggers that have any other triggers 's_pre' seconds before
        this_trig=remove_trig(all_trig_cell{i}, all_trig, s_pre);
        % find first and last frames for every trigger
        n_trig(i)=length(this_trig);
        this_id_start=nan(1,n_trig(i));
        this_id_end=this_id_start;
        for j=1:n_trig(i)
            % first frame (trig - s_pre)
            [~, this_id_start(j)] = min(abs(time_im - (this_trig(j) - s_pre)));
            this_id_start(j)=this_id_start(j)+1;
            % last frame (trig + s_post or 1 imaging frame before next trig);
            next_trig_id=find(all_trig-this_trig(j)>0,1);
            if isempty(next_trig_id)
                next_trig=inf;
            else
                next_trig=all_trig(next_trig_id)-im_dt;
            end
            this_id_end(j)=min(next_trig,this_trig(j) + s_post);
            [~, this_id_end(j)] = min(abs(time_im - this_id_end(j)));
        end
        % remove triggers that are too close to the next trigger
        tf_remove=time_im(this_id_end)<=this_trig;
        this_id_start(tf_remove)=[];
        this_id_end(tf_remove)=[];
        id_start{i} = this_id_start;
        id_end{i} = this_id_end;
        trig_length{i}=this_id_end-this_id_start+1;
        n_trig(i) = length(this_id_start);
    end  
        
    %% analysis for criteria (barcodes)
    % compute triggered avereages
    traces_bout_trig_crit = dm_compute_triggered_traces(bout1_starts,s_pre,s_post,dt,time_im,traces_interp);
    trace_motor_regr_bout_trig_crit = dm_compute_triggered_traces(bout1_starts,s_pre,s_post,dt,time_im,trace_motor_regr_interp);
    
    % compute scores
    scores_im = compute_scores_here(traces_bout_trig_crit,time_trig,clipping);
    scores_motor_regr = compute_scores_here(trace_motor_regr_bout_trig_crit,time_trig,clipping)';
    scores_be=trials.bout_duration.first;
    
    % compute criteria
    crit_im = nan(n_ROIs,4);
    temp_block2 = nanmean(scores_im(:,11:20),2);
    temp_block3 = nanmean(scores_im(:,21:30),2);
    temp_block7 = nanmean(scores_im(:,61:70),2);
    temp_block8 = nanmean(scores_im(:,71:80),2);
    temp_block12 = nanmean(scores_im(:,111:120),2);
    crit_im(:,1) = temp_block3 - temp_block2;
    crit_im(:,2) = temp_block7 - temp_block3;
    crit_im(:,3) = temp_block8 - temp_block7;
    crit_im(:,4) = temp_block12 - temp_block8;
    
    crit_motor_regr = nan(1,4);
    temp_block2 = nanmean(scores_motor_regr(11:20),2);
    temp_block3 = nanmean(scores_motor_regr(21:30),2);
    temp_block7 = nanmean(scores_motor_regr(61:70),2);
    temp_block8 = nanmean(scores_motor_regr(71:80),2);
    temp_block12 = nanmean(scores_motor_regr(111:120),2);
    crit_motor_regr(1) = temp_block3 - temp_block2;
    crit_motor_regr(2) = temp_block7 - temp_block3;
    crit_motor_regr(3) = temp_block8 - temp_block7;
    crit_motor_regr(4) = temp_block12 - temp_block8;

    crit_be=nan(1,4);
    temp_block2 = nanmean(scores_be(11:20));
    temp_block3 = nanmean(scores_be(21:30));
    temp_block7 = nanmean(scores_be(61:70));
    temp_block8 = nanmean(scores_be(71:80));
    temp_block12 = nanmean(scores_be(111:120));
    crit_be(1) = temp_block3 - temp_block2;
    crit_be(2) = temp_block7 - temp_block3;
    crit_be(3) = temp_block8 - temp_block7;
    crit_be(4) = temp_block12 - temp_block8;
  
    %% analysis for sensorimotor clustering
    % compute trigered averages
    scores=nan(n_ROIs,num_trig_types);
    mean_traces=nan(n_ROIs,n_frames_trig,num_trig_types,'single');
    ste_traces=mean_traces;
    for i=1:num_trig_types
        [mean_traces(:,:,i), ste_traces(:,:,i), scores(:,i)] = compute_trigaver_here(...
            id_start{i},...
            id_end{i},...
            trig_length{i},...
            n_ROIs,...
            n_frames_trig,...
            n_trig(i),...
            traces_interp,...
            baseline_tf,...
            score_tf(i,:));
    end
    traces_gr_trig_mean = mean_traces(:,:,1);
    traces_bout_trig_mean = mean_traces(:,:,2);
    traces_gr_trig_ste = ste_traces(:,:,1);
    traces_bout_trig_ste = ste_traces(:,:,2);
    scores_gr = scores(:,1);
    scores_bouts = scores(:,2);
    clear mean_traces ste_traces;
    
    %% bootstrapping
    % to build null ditributions for sens and motor scores
    % and for criteria
    chunk_start_ids = [1 : chunk_length : n_frames n_frames+1];
    n_chunks=length(chunk_start_ids(chunk_start_ids~=n_frames & chunk_start_ids~=n_frames+1));
    H0_scores=nan(n_ROIs,num_trig_types,n_boots);
    H0_im=nan(n_ROIs,n_boots);
    H0_motor_regr=nan(1,n_boots);
    H0_be=nan(1,n_boots);
    for ii=1:n_boots
        % for sensorymotor clustering
        chunk_ord=randperm(n_chunks);
        t=1;
        shuffled_traces=nan(n_ROIs,n_frames,'single');
        for j=1:n_chunks
            this_chunk_ids = chunk_start_ids(chunk_ord(j)) : chunk_start_ids(chunk_ord(j)+1)-1;
            this_chink_length = length(this_chunk_ids);
            shuffled_traces(:,t:t+this_chink_length-1) = traces_interp(:,this_chunk_ids);
            t=t+this_chink_length;
        end
        for i=1:num_trig_types
            H0_scores(:,i,ii) = compute_trigaver_here_fast(...
                id_start{i},...
                id_end{i},...
                trig_length{i},...
                n_ROIs,...
                n_frames_trig,...
                n_trig(i),...
                shuffled_traces,...
                baseline_tf,...
                score_tf(i,:));
        end
        % for criteria
        rand_ind = [1:10 randperm(110)+10];       
        scores_im_shuf=scores_im(:,rand_ind);
        temp_block2 = nanmean(scores_im_shuf(:,11:20),2);
        temp_block3 = nanmean(scores_im_shuf(:,21:30),2);
        H0_im(:,ii)=temp_block3-temp_block2;
        scores_motor_regr_shuf=scores_motor_regr(rand_ind);
        temp_block2 = nanmean(scores_motor_regr_shuf(11:20));
        temp_block3 = nanmean(scores_motor_regr_shuf(21:30));
        H0_motor_regr(ii)=squeeze(temp_block3-temp_block2);
        scores_be_shuf = scores_be(rand_ind);
        temp_block2 = nanmean(scores_be_shuf(11:20));
        temp_block3 = nanmean(scores_be_shuf(21:30));
        H0_be(ii)=temp_block3-temp_block2;
        progressbar([],[],ii/n_boots,[]);
    end

    %% find significant criteria
    p_thresh = 0.05;
    prcnt_min_im_5 = prctile(H0_im,p_thresh*100/2,2);
    prcnt_max_im_5 = prctile(H0_im,100-p_thresh*100/2,2);
    prcnt_min_motor_regr = prctile(H0_motor_regr,p_thresh*100/2);
    prcnt_max_motor_regr = prctile(H0_motor_regr,100-p_thresh*100/2);
    prcnt_min_be = prctile(H0_be,p_thresh*100/2);
    prcnt_max_be = prctile(H0_be,100-p_thresh*100/2);
    p_thresh = 0.15;
    prcnt_min_im_15 = prctile(H0_im,p_thresh*100/2,2);
    prcnt_max_im_15 = prctile(H0_im,100-p_thresh*100/2,2);

    crit_im_signif_5 = zeros(n_ROIs,4,'single');
    crit_im_signif_15 = zeros(n_ROIs,4,'single');
    for i=1:4
        crit_im_signif_5(crit_im(:,i)<prcnt_min_im_5,i)=-1;
        crit_im_signif_5(crit_im(:,i)>prcnt_max_im_5,i)=1;
        crit_im_signif_15(crit_im(:,i)<prcnt_min_im_15,i)=-1;
        crit_im_signif_15(crit_im(:,i)>prcnt_max_im_15,i)=1;
    end
    crit_motor_regr_signif = zeros(1,4,'single');
    crit_motor_regr_signif(crit_motor_regr<prcnt_min_motor_regr)=-1;
    crit_motor_regr_signif(crit_motor_regr>prcnt_max_motor_regr)=1;
    crit_be_signif = zeros(1,4);
    crit_be_signif(crit_be<prcnt_min_be)=-1;
    crit_be_signif(crit_be>prcnt_max_be)=1;
    clear H0_im H0_motor_regr H0_be
    
    %% cluster into sensory and motor cells
    thresh95=prctile(H0_scores,100-5,3);
    tf_sens=scores_gr>thresh95(:,1);
    tf_mot=scores_bouts>thresh95(:,2) & scores_gr<thresh95(:,1);
    sensmot_clust=zeros(size(scores,1),1,'uint8');
    sensmot_clust(tf_sens)=1;
    sensmot_clust(tf_mot)=2;
    clear tf_sens tf_mot;
     
    %% fit time constants to sensory ROIs   
    % build regressors for different time constants
    % convert times into indices
    gr_starts_ids=round(gr_starts/dt);
    gr_ends_ids=round(gr_ends/dt);
    % model activity of a neuron which is active during trials of forward
    % moving grating and silent otherwise
    model_act = zeros(1,n_frames);
    for i=1:length(gr_starts_ids)
        model_act(gr_starts_ids(i):gr_ends_ids(i))=1;
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
    traces_leaky_integrator = zeros(n_taus,n_frames);
    % if tau == dt, no integration
    for i = 1:n_taus
        tau=probing_taus(i);
        tau=dt/tau;
        this_trace_leaky_integrator = zeros(1,n_frames);
        if tau<inf
            for t=2:n_frames
                this_trace_leaky_integrator(t)=tau*model_act(t)-(tau-1)*this_trace_leaky_integrator(t-1);
            end
        else
            this_trace_leaky_integrator = model_act;
        end
        this_trace_leaky_integrator=conv(this_trace_leaky_integrator,my_kernel);
        traces_leaky_integrator(i,:) = this_trace_leaky_integrator(1:n_frames);
    end
    
    % find only sensory traces
    ids_sens = find(sensmot_clust==1);
    traces_interp_sens = traces_interp(ids_sens,:);
    n_ROIs_sens = length(ids_sens);

    % compute correlation values between real and modelled traces, find the best tau for each sensory ROI
    % do this for blocks of 10 trials
    corrcoeffs  = nan(n_ROIs_sens,n_taus,12);
    tt=0;
    for t=1:10:length(gr_starts_ids)
        tt=tt+1;
        ids = find(time_im>=gr_starts(t)-7.5 & time_im<gr_ends(t+9)+7.5);
        traces_interp_mean = nanmean(traces_interp_sens(:,ids),2);
        X = traces_interp_sens(:,ids) - traces_interp_mean;
        X_sqr = sum(X.^2,2,'omitnan');

        traces_leaky_integrator_mean = mean(traces_leaky_integrator(:,ids),2);
        Y = traces_leaky_integrator(:,ids) - traces_leaky_integrator_mean;
        Y_sqr = sum(Y.^2,2);
        
        for i = 1:n_taus
            corrcoeffs(:,i,tt) = sum(X.*Y(i,:),2,'omitnan')./...
                sqrt(X_sqr*Y_sqr(i));
        end
        progressbar([],[],[],tt/12);
    end
    [~,best_corrcoeff_ids]= max(corrcoeffs,[],2);
    time_constants_sens = squeeze(probing_taus(best_corrcoeff_ids));
    time_constants = nan(n_ROIs,12);
    time_constants(ids_sens,:) = time_constants_sens;
    clear traces_leaky_integrator_mean model_act best_corrcoeff_ids traces_leaky_integrator this_trace_leaky_integrator traces_interp_sens traces_interp_mean;

    %% save
    save(filename,...
    'time_trig', ...
    'traces_bout_trig_crit', ...
    'trace_motor_regr_bout_trig_crit', ...
    'traces_gr_trig_mean', ...
    'traces_bout_trig_mean', ...
    'traces_gr_trig_ste', ...
    'traces_bout_trig_ste', ...
    'sensmot_clust', ...
    'scores_be', ...
    'scores_im', ...
    'scores_motor_regr', ...
    'crit_be', ...
    'crit_be_signif', ...
    'crit_im', ...
    'crit_im_signif_5', ...
    'crit_im_signif_15', ...
    'crit_motor_regr', ...
    'crit_motor_regr_signif', ...
    'time_constants' ...
    );
    
progressbar(f/n_fish,[],[],[]);
end


function [this_trig] = remove_trig(this_trig, all_trig, s_pre)
num_trig=length(this_trig);
remove_trig_tf = false(1,num_trig);
for i=1:num_trig
    all_dist=all_trig-this_trig(i);
    if any(all_dist<0 & -all_dist<s_pre)
        remove_trig_tf(i)=true;
    end
end
this_trig(remove_trig_tf)=[];
end

function [mean_traces, ste_traces, scores] = compute_trigaver_here(...
    this_id_start,...
    this_id_end,...
    this_trig_length,...
    n_ROIs,...
    n_frames,...
    n_trig,...
    traces,...
    baseline_tf,...
    score_tf)
trig_traces=nan(n_ROIs,n_frames,n_trig,'single');
for j=1:n_trig
    trig_traces(:,1:this_trig_length(j),j) = traces(:,this_id_start(j):this_id_end(j));
end
trig_traces=trig_traces-nanmean(trig_traces(:,baseline_tf,:),2);
mean_traces=nanmean(trig_traces,3);
ste_traces=nanstd(trig_traces,[],3)./sqrt(sum(~isnan(trig_traces),3));
scores=nanmean(mean_traces(:,score_tf),2);
end

function scores = compute_trigaver_here_fast(...
    this_id_start,...
    this_id_end,...
    this_trig_length,...
    n_ROIs,...
    n_frames,...
    n_trig,...
    traces,...
    baseline_tf,...
    score_tf)
trig_traces=nan(n_ROIs,n_frames,n_trig,'single');
for j=1:n_trig
    trig_traces(:,1:this_trig_length(j),j) = traces(:,this_id_start(j):this_id_end(j));
end
trig_traces=trig_traces-nanmean(trig_traces(:,baseline_tf,:),2);
mean_traces=nanmean(trig_traces,3);
scores=nanmean(mean_traces(:,score_tf),2);
end

function [scores] = compute_scores_here(data,time,clipping)
scores=squeeze(nanmean(data(:,time>0 & time<=clipping,:),2));
end