% this progam analyses behavior of the whole-brain imaging experiment.
% It perfoms the same computations as behavioral_analysis.m
% but using a different sourse data format

%% pathnames
clc; close all; clear all;
pathname_raw_data='J:\_Shared\experiments\E0031_acute_adaptation\v02_lightsheet\data_behavior\';
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_behavior=[pathname_imaging_results 'behavior\'];

%% some important variables
dt=0.005;
spont_dur=120;
rest_dur=7.5;
trial_dur=15;
full_trial_dur=trial_dur+2*rest_dur;
calib_notr=0; %10;
pre_notr=10;
main_notr=40;
post_notr=0;

spont_notr=spont_dur/full_trial_dur;
total_notr=calib_notr+pre_notr+main_notr+post_notr;
t_pre_start=full_trial_dur*calib_notr+spont_dur;
t_main_start=full_trial_dur*(calib_notr+pre_notr)+spont_dur;
t_post_start=full_trial_dur*(calib_notr+pre_notr+main_notr)+spont_dur;
gr_starts=spont_dur+(rest_dur:full_trial_dur:rest_dur+(total_notr-1)*full_trial_dur);
gr_ends=gr_starts+trial_dur;
time_be=dt:dt:total_notr*30+spont_dur;

bout_par_names={'bout_duration','next_interbout_duration','sum_bout_vigor','mean_bout_freq','median_bout_freq','max_bout_freq','mean_bout_amp','median_bout_amp','max_bout_amp','distance_bout','distance_bout_interbout'};
trial_par_names={'number_of_bouts','latency'};
bout_arrays_names={'vigor','power','grspeed','flick_amps','flick_freqs'};
measures={'first','mean','median','max','std','ste'};


%% loop through all fish
cd(pathname_raw_data);
all_fish=dm_dir('*_f*');
num_fish=length(all_fish);
for f=1:num_fish
    this_fish=all_fish{f};
    filename_data=[this_fish '_data'];
    %% metadata
    path_to_fish=[pathname_raw_data this_name '\'];
    cd(path_to_fish);
    meta=load(filename_data,'meta');
    meta=meta.meta;
    time_exp_dur=datetime(meta.general.t_protocol_end([1:10 12:19]),'InputFormat','yyyy-MM-ddHH:mm:ss')-datetime(meta.general.t_protocol_start([1:10 12:19]),'InputFormat','yyyy-MM-ddHH:mm:ss');
    comments=meta.general.animal.comments;
    
    %% stimulus
    % read the values
    stim=load(filename_data,'stim');
    stim=stim.stim;
    if isfield(stim, 'grspeed_premainpost')
        stim=rmfield(stim, {'grspeed_premainpost','base_grspeed_premainpost','swim_premainpost'});
    end
    
    bad_points=diff(stim.time)==0;
    if sum(bad_points)>0
        for i=fieldnames(stim)'
            stim.(cell2mat(i))(bad_points)=[];
        end
    end
    tf_prepost=(stim.time>=t_pre_start & stim.time<t_main_start) | stim.time>t_post_start;
    tf_main=stim.time>=t_main_start & stim.time<t_post_start;
    grspeed=zeros(1,length(stim.time));
    grspeed(tf_prepost)=stim.grspeed_prepost(tf_prepost);
    grspeed(tf_main)=stim.grspeed_main(tf_main);
    base_grspeed=zeros(1,length(stim.time));
    base_grspeed(tf_prepost)=stim.base_grspeed_prepost(tf_prepost);
    base_grspeed(tf_main)=stim.base_grspeed_main(tf_main);
    swim=zeros(1,length(stim.time));
    swim(tf_prepost)=stim.swim_prepost(tf_prepost);
    swim(tf_main)=stim.swim_main(tf_main);
    
    % find bout starts and ends (detected in real time during the experiment)
    bout_starts_ids=find(diff([0 swim])==1);
    bout_ends_ids=find(diff([0 swim 0])==-1);
    if ~isempty(bout_ends_ids)
        if bout_ends_ids(end)>length(stim.time)
            bout_ends_ids(end)=length(stim.time);
        end
    end
    nob=length(bout_starts_ids);
    % find reafference of all these bouts
    bouts=[];
    gain=ones(1,length(stim.time));
    gain(tf_main)=stim.gain_main(tf_main);
    lag=zeros(1,length(stim.time));
    lag(tf_main)=stim.lag_main(tf_main);
    gd_starts=zeros(1,length(stim.time));
    gd_starts(tf_main)=stim.gd_starts_main(tf_main);
    gd_ends=zeros(1,length(stim.time));
    gd_ends(tf_main)=stim.gd_ends_main(tf_main);
    bouts.gain=gain(bout_starts_ids);
    bouts.lag=lag(bout_starts_ids);
    bouts.measured_lag=nan(1,nob);
    bouts.gain_drop_start=gd_starts(bout_starts_ids);
    bouts.gain_drop_start(isnan(bouts.gain_drop_start))=0;
    bouts.gain_drop_end=gd_ends(bout_starts_ids);
    bouts.gain_drop_end(isnan(bouts.gain_drop_end))=0;
    bouts.gain_drop_duration=round((bouts.gain_drop_end-bouts.gain_drop_start)*1000)/1000;
    notr=sum(diff(base_grspeed)==-10 & base_grspeed(1:end-1)==10);
    
    % prepare conditions
    cond = find_all_conditions (bouts);  
    
    %% behavioral traces
    % read the values
    be=load(filename_data,'be');
    be=be.be;
    % correct the repeated time points
    bad_points=diff(be.time)==0;
    if sum(bad_points)>0
        be.time(bad_points)=[];
        be.tail(bad_points)=[];
    end
    % interpolate tail trace, swim and grating speed
    tail=interp1(be.time,be.tail,time_be);
    swim=round(interp1(stim.time,swim,time_be))==1;
    grspeed=interp1(stim.time,grspeed,time_be);
    % scale the tail
    tail=tail-nanmean(tail(~swim));
    tail=tail/nanstd(tail(swim));
    tail(isnan(tail))=0;
    % build vigor trace
    vigor=nan(1,length(tail));
    for i=0.05/dt:length(tail)
        vigor(i)=std(tail(i-(0.05/dt-1):i));
    end
    if length(vigor)<total_notr*full_trial_dur/dt % (for some fish the very last frame is not there)
        vigor=[vigor nan(1,total_notr*full_trial_dur/dt-length(vigor))];
    end
    
    %% detect bouts accurately (and find min and max tail positions)
    bouts.start=nan(1,nob);
    bouts.end=nan(1,nob);
    bouts.bad_bouts=false(1,nob);
    mm_time_array=nan(nob,1000);
    mm_val_array=mm_time_array;
    max_tf_array=false(nob,1000);
    for b=1:nob
        bst0=round((stim.time(bout_starts_ids(b))-0.1)/dt);
        bet0=round(stim.time(bout_ends_ids(b))/dt);
        bouts.start(b)=stim.time(bout_starts_ids(b));
        bouts.end(b)=stim.time(bout_ends_ids(b));
        
        % find individual flicks
        this_mm_time_array=[];
        this_mm_val_array=[];
        this_max_tf_array=[];
        c=0;
        if bst0<2
            bst0=2;
        end
        if bet0>total_notr*full_trial_dur/dt-1
            bet0=total_notr*full_trial_dur/dt-1;
        end
        for i=bst0:bet0
            if tail(i)>=tail(i-1) && tail(i)>tail(i+1)
                c=c+1;
                this_mm_time_array(c)=time_be(i);
                this_mm_val_array(c)=tail(i);
                this_max_tf_array(c)=true;
            elseif tail(i)<=tail(i-1) && tail(i)<tail(i+1)
                c=c+1;
                this_mm_time_array(c)=time_be(i);
                this_mm_val_array(c)=tail(i);
                this_max_tf_array(c)=false;
            end
        end
        % remove all bad flicks
        if length(this_mm_time_array)>=2
            bad_flicks=[abs(diff(this_mm_val_array))<0.14 false] & [false flip(abs(diff(flip(this_mm_val_array)))<0.14)];
            if abs(this_mm_val_array(2)-this_mm_val_array(1))<0.14
                bad_flicks(1)=true;
            end
            if abs(this_mm_val_array(end)-this_mm_val_array(end-1))<0.14
                bad_flicks(end)=true;
            end
            bad_flicks=bad_flicks | [false flip(abs(diff(flip(this_mm_time_array)))>0.1)];
            this_mm_time_array(bad_flicks)=[];
            this_mm_val_array(bad_flicks)=[];
            this_max_tf_array(bad_flicks)=[];
            bad_flicks=diff(this_max_tf_array)==0;
            this_mm_time_array(bad_flicks)=[];
            this_mm_val_array(bad_flicks)=[];
            this_max_tf_array(bad_flicks)=[];
            if length(this_mm_time_array)>=2
                bouts.start(b)=this_mm_time_array(1);
                bouts.end(b)=this_mm_time_array(end);
            else
                bouts.bad_bouts(b)=true;
            end
        else
            bouts.bad_bouts(b)=true;
        end
        n=length(this_mm_time_array);
        if n>1000
            n=1000;
        end
        mm_time_array(b,1:n)=this_mm_time_array(1:n);
        mm_val_array(b,1:n)=this_mm_val_array(1:n);
        max_tf_array(b,1:n)=this_max_tf_array(1:n);
    end
    % shorten the flicks arrays
    mm_n=find(all(isnan(mm_time_array),1),1)-1;
    if ~isempty(mm_n)
        if mm_n>=50
            mm_n=50;
        end
        mm_time_array=mm_time_array(:,1:mm_n);
        mm_val_array=mm_val_array(:,1:mm_n);
        max_tf_array=max_tf_array(:,1:mm_n);
    end
    % create array of bad bouts
    bouts.long_bouts=bouts.end-bouts.start>=0.3;
    bouts.short_bouts=bouts.end-bouts.start<0.1;
    bouts.bad_bouts=bouts.bad_bouts | ... % bouts that are already bad (failed to detect them properly)
        bouts.short_bouts |... % bouts which are shorter than 100 ms
        [false bouts.start(2:end)-bouts.end(1:end-1)<0.1] | [bouts.start(2:end)-bouts.end(1:end-1)<0.1 false] |... % bouts with interbouts shorter than 100 ms (tipically, these are bouts which were detected as 2 bouts)
        all(isnan(mm_time_array),2)' |...  % bouts with no flicks
        max(diff(mm_time_array,1,2),[],2)'>0.1; % bouts with max delta flick time > 100 ms (this happens during some wierd ugly bouts)
    % find spontaneous bouts (i.e. bouts which started after trial start and finished before trial end)
    bouts.spont_bouts=true(1,nob);
    for g=1:notr
        bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
        bouts.spont_bouts(bouts_in_this_trial)=false;
    end
    % find bouts during the main part of the experiment (adaptation phase)
    bouts.main_bouts=stim.time(bout_ends_ids)<t_post_start & stim.time(bout_starts_ids)>=t_main_start;
    bouts.main_bouts=bouts.main_bouts';
    % good bouts: i.e. they are not bad, not spontaneous and happened during main part
    bouts.good_bouts=bouts.main_bouts & ~bouts.bad_bouts & ~bouts.spont_bouts;
    
    %% find bout parameters
    bouts.trial=nan(1,nob);
    bouts.bout_duration=bouts.end-bouts.start;
    bouts.next_interbout_duration=nan(1,nob);
    bouts.distance_bout=nan(1,nob);
    bouts.distance_bout_interbout=nan(1,nob);
    for g=1:notr
        bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
        for b=bouts_in_this_trial
            bouts.trial(b)=g;
            bs=round((bouts.start(b))/dt);
            be=round((bouts.end(b))/dt);
            bouts.distance_bout(b)=sum(grspeed(bs:be))*dt;
            if b<bouts_in_this_trial(end)
                bouts.next_interbout_duration(b)=bouts.start(b+1)-bouts.end(b);
                be=round((bouts.start(b+1))/dt);
                bouts.distance_bout_interbout(b)=sum(grspeed(bs:be))*dt;
            end
            bouts.sum_bout_vigor(b)=sum(vigor(round(bouts.start(b)/dt):round(bouts.end(b)/dt)))*dt;
            freq=1./(diff(mm_time_array(b,:))*2);
            bouts.mean_bout_freq(b)=nanmean(freq);
            bouts.median_bout_freq(b)=nanmedian(freq);
            bouts.max_bout_freq(b)=max(freq);
            amp=abs(diff(mm_val_array(b,:)));
            bouts.mean_bout_amp(b)=nanmean(amp);
            bouts.median_bout_amp(b)=nanmedian(amp);
            bouts.max_bout_amp(b)=max(amp);
            bs=round((bouts.start(b)-0.1)/dt);
            be=round(bouts.end(b)/dt);
            if be-bs>=1.1/dt
                be=bs+1.1/dt-1;
            end
            bouts.tail(b,:)=zeros(1,1.1/dt); % changed from nan
            bouts.tail(b,1:be-bs+1)=tail(bs:be);
            bouts.vigor(b,:)=zeros(1,1.1/dt); % changed from nan
            bouts.vigor(b,1:be-bs+1)=vigor(bs:be);
            bouts.power(b,:)=zeros(1,1.1/dt); % changed from nan
            bouts.power(b,1:be-bs+1)=tail(bs:be).^2;
            bouts.grspeed(b,:)=nan(1,1.1/dt);
            bouts.grspeed(b,1:be-bs+1)=grspeed(bs:be);
            bouts.flick_times(b,:)=mm_time_array(b,:);
            bouts.flick_values(b,:)=mm_val_array(b,:);
            bouts.flick_amps(b,:)=abs(diff(bouts.flick_values(b,:),1,2));
            bouts.flick_freqs(b,:)=1./(diff(bouts.flick_times(b,:),1,2)*2);
        end
    end
    for p=1:length(bout_par_names)
        bouts.(bout_par_names{p})(bouts.bad_bouts | bouts.spont_bouts)=nan;
    end
    for p=1:length(bout_arrays_names)
        bouts.(bout_arrays_names{p})(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    end
    bouts.tail(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    bouts.flick_times(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    bouts.flick_values(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    
    %% find trial averages
    trials=[];
    vigor_for_imagesc=nan(total_notr+spont_notr,full_trial_dur/dt);
    vigor_for_imagesc(1:spont_notr,:)=reshape(vigor(1:spont_dur/dt),full_trial_dur/dt,spont_notr)';
    for g=1:notr
        trial_start=round((gr_starts(g)-rest_dur)/dt)+1;
        trial_end=round((gr_ends(g)+rest_dur)/dt);
        vigor_for_imagesc(g+spont_notr,:)=vigor(trial_start:trial_end);
        bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
        trials.bouts_in_this_trial(g)={bouts_in_this_trial};
        trials.number_of_bouts(g)=length(bouts_in_this_trial);
        if ~isempty(bouts_in_this_trial)
            trials.latency(g)=bouts.start(bouts_in_this_trial(1))-gr_starts(g);
        else
            trials.latency(g)=nan;
        end
        for p=1:length(bout_par_names)
            this_par=bout_par_names{p};
            these_bouts=bouts.(this_par)(bouts_in_this_trial);
            these_bouts=these_bouts(~isnan(these_bouts));
            trials.(this_par).n(g)=length(these_bouts);
            if trials.(this_par).n(g)>0
                trials.(this_par).first(g)=these_bouts(1);
                trials.(this_par).mean(g)=mean(these_bouts);
                trials.(this_par).median(g)=median(these_bouts);
                trials.(this_par).max(g)=max(these_bouts);
                trials.(this_par).std(g)=std(these_bouts);
                trials.(this_par).ste(g)=trials.(this_par).std(g)/sqrt(n);
            else
                trials.(this_par).first(g)=nan;
                trials.(this_par).mean(g)=nan;
                trials.(this_par).median(g)=nan;
                trials.(this_par).max(g)=nan;
                trials.(this_par).std(g)=nan;
                trials.(this_par).ste(g)=nan;
            end
        end
        for p=1:length(bout_arrays_names)
            this_par=bout_arrays_names{p};
            these_bouts=bouts.(this_par)(bouts_in_this_trial,:);
            these_bouts=these_bouts(~all(isnan(these_bouts),2),:);
            trials.(this_par).n(g)=size(these_bouts,1);
            if trials.(this_par).n(g)>0
                trials.(this_par).first(g,:)=these_bouts(1,:);
                trials.(this_par).mean(g,:)=nanmean(these_bouts,1);
                trials.(this_par).median(g,:)=nanmedian(these_bouts,1);
                trials.(this_par).max(g,:)=nanmax(these_bouts,[],1);
                trials.(this_par).std(g,:)=nanstd(these_bouts,0,1);
                trials.(this_par).ste(g,:)=trials.(this_par).std(g,:)/sqrt(n);
            else
                trials.(this_par).first(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).mean(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).median(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).max(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).std(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).ste(g,:)=nan(1,size(bouts.(this_par),2));
            end
        end
    end
    
    %% find phase mean
    for g=1:floor(notr/10)
        these_trials=(g-1)*10+1:g*10;
        for p=1:length(trial_par_names)
            this_par=trial_par_names{p};
            these_bouts=trials.(this_par)(these_trials);
            these_bouts=these_bouts(~isnan(these_bouts));
            trials.phase_mean.(this_par).n(g)=length(these_bouts);
            if trials.phase_mean.(this_par).n(g)>0
                trials.phase_mean.(this_par).mean(g)=mean(these_bouts);
                trials.phase_mean.(this_par).median(g)=median(these_bouts);
                trials.phase_mean.(this_par).max(g)=max(these_bouts);
                trials.phase_mean.(this_par).std(g)=std(these_bouts);
                trials.phase_mean.(this_par).ste(g)=trials.phase_mean.(this_par).std(g)/sqrt(trials.phase_mean.(this_par).n(g));
            else
                trials.phase_mean.(this_par).mean(g)=nan;
                trials.phase_mean.(this_par).median(g)=nan;
                trials.phase_mean.(this_par).max(g)=nan;
                trials.phase_mean.(this_par).std(g)=nan;
                trials.phase_mean.(this_par).ste(g)=nan;
            end
        end
        for p=1:length(bout_par_names)
            this_par=bout_par_names{p};
            for m=1:length(measures)
                this_measure_name=measures{m};
                these_bouts=trials.(this_par).(this_measure_name)(these_trials);
                these_bouts=these_bouts(~isnan(these_bouts));
                trials.phase_mean.(this_par).(this_measure_name).n(g)=length(these_bouts);
                if trials.phase_mean.(this_par).(this_measure_name).n(g)>0
                    trials.phase_mean.(this_par).(this_measure_name).mean(g)=mean(these_bouts);
                    trials.phase_mean.(this_par).(this_measure_name).median(g)=median(these_bouts);
                    trials.phase_mean.(this_par).(this_measure_name).max(g)=max(these_bouts);
                    trials.phase_mean.(this_par).(this_measure_name).std(g)=std(these_bouts);
                    trials.phase_mean.(this_par).(this_measure_name).ste(g)=trials.phase_mean.(this_par).(this_measure_name).std(g)/sqrt(trials.phase_mean.(this_par).(this_measure_name).n(g));
                else
                    trials.phase_mean.(this_par).(this_measure_name).mean(g)=nan;
                    trials.phase_mean.(this_par).(this_measure_name).median(g)=nan;
                    trials.phase_mean.(this_par).(this_measure_name).max(g)=nan;
                    trials.phase_mean.(this_par).(this_measure_name).std(g)=nan;
                    trials.phase_mean.(this_par).(this_measure_name).ste(g)=nan;
                end
            end
        end
        for p=1:length(bout_arrays_names)
            this_par=bout_arrays_names{p};
            for m=1:length(measures)
                this_measure_name=measures{m};
                these_bouts=trials.(this_par).(this_measure_name)(these_trials,:);
                these_bouts=these_bouts(~all(isnan(these_bouts),2),:);
                trials.phase_mean.(this_par).(this_measure_name).n(g)=size(these_bouts,1);
                if trials.phase_mean.(this_par).(this_measure_name).n(g)>0
                    trials.phase_mean.(this_par).(this_measure_name).mean(g,:)=nanmean(these_bouts,1);
                    trials.phase_mean.(this_par).(this_measure_name).median(g,:)=nanmedian(these_bouts,1);
                    trials.phase_mean.(this_par).(this_measure_name).max(g,:)=nanmax(these_bouts,[],1);
                    trials.phase_mean.(this_par).(this_measure_name).std(g,:)=nanstd(these_bouts,0,1);
                    trials.phase_mean.(this_par).(this_measure_name).ste(g,:)=trials.phase_mean.(this_par).(this_measure_name).std(g,:)/sqrt(trials.phase_mean.(this_par).(this_measure_name).n(g));
                else
                    trials.phase_mean.(this_par).(this_measure_name).mean(g,:)=nan(1,size(bouts.(this_par),2));
                    trials.phase_mean.(this_par).(this_measure_name).median(g,:)=nan(1,size(bouts.(this_par),2));
                    trials.phase_mean.(this_par).(this_measure_name).max(g,:)=nan(1,size(bouts.(this_par),2));
                    trials.phase_mean.(this_par).(this_measure_name).std(g,:)=nan(1,size(bouts.(this_par),2));
                    trials.phase_mean.(this_par).(this_measure_name).ste(g,:)=nan(1,size(bouts.(this_par),2));
                end
            end
        end
    end
    
    %% find condtition averages
    for c=1:length(cond.names)
        this_cond_name=cond.names{c};
        this_cond_vals=cond.vals{c};
        bouts.condition_mean.(this_cond_name).conditions=cond.xlabel{c};
        for i=1:length(bouts.condition_mean.(this_cond_name).conditions)
            tf=bouts.good_bouts & bouts.lag==this_cond_vals(1,i) & bouts.gain==this_cond_vals(2,i) & bouts.gain_drop_start==this_cond_vals(3,i) & bouts.gain_drop_end==this_cond_vals(4,i);
            for p=1:length(bout_par_names)
                this_par=bout_par_names{p};
                these_bouts=bouts.(this_par)(tf);
                these_bouts=these_bouts(~isnan(these_bouts));
                bouts.condition_mean.(this_cond_name).(this_par).n(i)=length(these_bouts);
                if bouts.condition_mean.(this_cond_name).(this_par).n(i)>0
                    bouts.condition_mean.(this_cond_name).(this_par).mean(i)=mean(these_bouts);
                    bouts.condition_mean.(this_cond_name).(this_par).median(i)=median(these_bouts);
                    bouts.condition_mean.(this_cond_name).(this_par).max(i)=max(these_bouts);
                    bouts.condition_mean.(this_cond_name).(this_par).std(i)=std(these_bouts);
                    bouts.condition_mean.(this_cond_name).(this_par).ste(i)=bouts.condition_mean.(this_cond_name).(this_par).std(i)/sqrt(bouts.condition_mean.(this_cond_name).(this_par).n(i));
                else
                    bouts.condition_mean.(this_cond_name).(this_par).mean(i)=nan;
                    bouts.condition_mean.(this_cond_name).(this_par).median(i)=nan;
                    bouts.condition_mean.(this_cond_name).(this_par).max(i)=nan;
                    bouts.condition_mean.(this_cond_name).(this_par).std(i)=nan;
                    bouts.condition_mean.(this_cond_name).(this_par).ste(i)=nan;
                end
            end
            for p=1:length(bout_arrays_names)
                this_par=bout_arrays_names{p};
                these_bouts=bouts.(this_par)(tf,:);
                these_bouts=these_bouts(~all(isnan(these_bouts),2),:);
                bouts.condition_mean.(this_cond_name).(this_par).n(i)=size(these_bouts,1);
                if bouts.condition_mean.(this_cond_name).(this_par).n(i)
                    bouts.condition_mean.(this_cond_name).(this_par).mean(i,:)=nanmean(these_bouts,1);
                    bouts.condition_mean.(this_cond_name).(this_par).median(i,:)=nanmedian(these_bouts,1);
                    bouts.condition_mean.(this_cond_name).(this_par).max(i,:)=nanmax(these_bouts,[],1);
                    bouts.condition_mean.(this_cond_name).(this_par).std(i,:)=nanstd(these_bouts,0,1);
                    bouts.condition_mean.(this_cond_name).(this_par).ste(i,:)=bouts.condition_mean.(this_cond_name).(this_par).std(i,:)/sqrt(bouts.condition_mean.(this_cond_name).(this_par).n(i));
                else
                    bouts.condition_mean.(this_cond_name).(this_par).mean(i,:)=nan(1,size(bouts.(this_par),2));
                    bouts.condition_mean.(this_cond_name).(this_par).median(i,:)=nan(1,size(bouts.(this_par),2));
                    bouts.condition_mean.(this_cond_name).(this_par).max(i,:)=nan(1,size(bouts.(this_par),2));
                    bouts.condition_mean.(this_cond_name).(this_par).std(i,:)=nan(1,size(bouts.(this_par),2));
                    bouts.condition_mean.(this_cond_name).(this_par).ste(i,:)=nan(1,size(bouts.(this_par),2));
                end
            end
        end
    end
    
    %% save the data
    filename=[pathname_behavior this_fish '_behavior.mat'];
    save(filename,'time_be','tail','grspeed','swim','bouts','trials','meta');
end