%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_behavior=[pathname_imaging_results 'behavior\'];
pathname_traces=[pathname_imaging_results 'traces\'];
pathname_trigaver=[pathname_imaging_results 'triggered_traces\'];
pathname_clustering=[pathname_imaging_results 'clustering\'];

%% timing parameters
time_pre=1;
time_post=4;
dt=0.05;
n_frames=(time_pre+time_post)/dt;
time_trig=round((-time_pre+dt:dt:time_post)*100)/100;
score_tf=[time_trig>0;time_trig>0 & time_trig<2]; % sens tf; mot tf
baseline_tf=1:time_pre/dt-1;

%% parameters for bootstrapping
n_bootstrap = 1000;
chunk_length = 23; % seconds
chunk_length = round(chunk_length/dt); % frames

%% loop through all fish
all_fish=dm_dir([pathname_traces '*_f*_traces.mat']);
all_fish=strrep(all_fish,'_traces.mat','');
n_fish=length(all_fish);
progressbar('Fish progress...','Interpolating traces...','Computing H0 for scores')
for f=4:n_fish
    fish_id=all_fish{f};
    filename_trigaver = [pathname_trigaver fish_id '_trig_traces.mat'];
    filename_clustering = [pathname_clustering fish_id '_clustering.mat'];
    
    %% load the data
    load([pathname_traces fish_id '_traces.mat'],'traces','time_offsets');
    n_ROIs=size(traces,1);
    load([pathname_behavior fish_id '_behavior.mat'],'bouts','meta');
    im_dt=1/meta.fs;
    time_im_ori=im_dt:im_dt:size(traces,2)*im_dt;
    time_im=round((-time_pre+dt:dt:time_im_ori(end)+time_post)*100)/100;
    num_frames_total=length(time_im);
    
    %% create triggers
    % sensory triggers
    forw_gr_starts = 120+7.5:30:time_im_ori(end);
    forw_gr_ends=forw_gr_starts+15;
    rev_gr_starts = 120+2.5:30:time_im_ori(end);
    % extract motor triggeres
    bout_starts=bouts.start;
    % build all triggers together
    all_trig=sort([forw_gr_starts, forw_gr_ends, rev_gr_starts, bout_starts]);
    all_trig_cell = {forw_gr_starts, bout_starts};
    num_trig_types=length(all_trig_cell);
    % work with triggers
    n_trig=zeros(1,num_trig_types);
    id_start=cell(1,num_trig_types);
    id_end=id_start;
    trig_length=id_start;
    for i=1:num_trig_types
        % remove triggers that have any other triggers 'time_pre' seconds before
        this_trig=remove_trig(all_trig_cell{i}, all_trig, time_pre);
        % find first and last frames for every trigger
        n_trig(i)=length(this_trig);
        this_id_start=nan(1,n_trig(i));
        this_id_end=this_id_start;
        for j=1:n_trig(i)
            % first frame (trig - time_pre)
            [~, this_id_start(j)] = min(abs(time_im - (this_trig(j) - time_pre)));
            this_id_start(j)=this_id_start(j)+1;
            % last frame (trig + time_post or 1 imaging frame before next trig);
            next_trig_id=find(all_trig-this_trig(j)>0,1);
            if isempty(next_trig_id)
                next_trig=inf;
            else
                next_trig=all_trig(next_trig_id)-im_dt;
            end
            this_id_end(j)=min(next_trig,this_trig(j) + time_post);
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
    
    %% interpolate all traces to time_im
    traces_interp = nan(n_ROIs,num_frames_total,'single');
    for i=1:n_ROIs
        traces_interp(i,:)=interp1(time_im_ori+time_offsets(i),traces(i,:),time_im);
        progressbar([],i/n_ROIs,[]);
    end
    
    %% compute triggered averages and sensory and motor scores
    scores=nan(n_ROIs,num_trig_types);
    mean_traces=nan(n_ROIs,n_frames,num_trig_types,'single');
    ste_traces=mean_traces;
    for i=1:num_trig_types
        [mean_traces(:,:,i), ste_traces(:,:,i), scores(:,i)] = compute_trigaver_here(...
            id_start{i},...
            id_end{i},...
            trig_length{i},...
            n_ROIs,...
            n_frames,...
            n_trig(i),...
            traces_interp,...
            baseline_tf,...
            score_tf(i,:));
    end
     
    %% build null ditributions for sens and motor scores
    chunk_start_ids = [1 : chunk_length : num_frames_total num_frames_total+1];
    n_chunks=length(chunk_start_ids(chunk_start_ids~=num_frames_total & chunk_start_ids~=num_frames_total+1));
    H0_scores=nan(n_ROIs,num_trig_types,n_bootstrap);
    for ii=1:n_bootstrap
        chunk_ord=randperm(n_chunks);
        t=1;
        shuffled_traces=nan(n_ROIs,num_frames_total,'single');
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
                n_frames,...
                n_trig(i),...
                shuffled_traces,...
                baseline_tf,...
                score_tf(i,:));
        end
        progressbar([],[],ii/n_bootstrap);
    end
    
    %% prepare relevant data for saving
    traces_gr_trig_mean = mean_traces(:,:,1);
    traces_bout_trig_mean = mean_traces(:,:,2);
    traces_gr_trig_ste = ste_traces(:,:,1);
    traces_bout_trig_ste = ste_traces(:,:,2);
    scores_gr = scores(:,1);
    scores_bouts = scores(:,2);
    
    %% cluster into sensory and motor cells
    thresh95=prctile(H0_scores,100-5,3);
    tf_sens=scores_gr>thresh95(:,1);
    tf_mot=scores_bouts>thresh95(:,2) & scores_gr<thresh95(:,1);
    sensmot_clust=zeros(size(scores,1),1,'uint8');
    sensmot_clust(tf_sens)=1;
    sensmot_clust(tf_mot)=2;
          
    %% save everything for this fish
    save(filename_trigaver,'traces_gr_trig_mean','traces_bout_trig_mean','traces_gr_trig_ste','traces_bout_trig_ste','time_trig');
    save(filename_clustering,'sensmot_clust');
    save(strrep(filename_trigaver,'_trig_traces.mat','_H0_delete_this.mat'),'H0_scores');
    
    progressbar(f/n_fish,[],[]);
end

function [this_trig] = remove_trig(this_trig, all_trig, time_pre)
num_trig=length(this_trig);
remove_trig_tf = false(1,num_trig);
for i=1:num_trig
    all_dist=all_trig-this_trig(i);
    if any(all_dist<0 & -all_dist<time_pre)
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