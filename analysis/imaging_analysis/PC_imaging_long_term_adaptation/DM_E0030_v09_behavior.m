% this progam analyses behavior of the long-term adaptation experiment during
% PC lightsheet imaging. It perfoms the same computations as behavioral_analysis.m
% but using a different sourse data format

%% pathnames
clc; close all; clear all;
% original path to the raw behavioral hdf5 data saved by Stytra
pathname_data='J:\_Shared\experiments\E0030_long_term_adaptation\v09_lightsheet\data_behavior\';
% final path to pre-processed data
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_behavior=[pathname_MS_data 'PC_imaging_long_term_adaptation\behavior\'];

%% timing variables
dt=0.005;
rest_dur=7.5;
trial_dur=15;
full_trial_dur=trial_dur+2*rest_dur;
calib_notr=10;
pre_notr=10;
main_notr=50;
post_notr=50;
notr=calib_notr+pre_notr+main_notr+post_notr;
exp_dur=notr*full_trial_dur;
gr_starts=rest_dur:full_trial_dur:rest_dur+(notr-1)*full_trial_dur;
gr_ends=gr_starts+trial_dur;
t_pre_start=full_trial_dur*calib_notr;
t_main_start=full_trial_dur*(calib_notr+pre_notr);
t_post_start=full_trial_dur*(calib_notr+pre_notr+main_notr);
time_be=dt:dt:exp_dur;

%% bout parameters
bout_par_names={'bout_duration','next_interbout_duration'};
trial_par_names={'number_of_bouts','latency'};
bout_arrays_names={'power','grspeed'};
measures={'first','mean'};

%% loop through all fish
cd(pathname_data);
all_fish=dm_dir('*_f*');
num_fish=length(all_fish);
for f=1:num_fish
    this_fish=all_fish{f};
    this_path=[pathname_data this_fish '\'];
    cd(this_path);
    filename=[pathname_behavior this_fish '_behavior.mat'];
    
    %% work with metadata
    filename_meta=dm_dir('*_metadata.json');
    meta2=read_json(filename_meta{end});
    meta=struct;
    meta.dpf=meta2.general.animal.age;
    meta.genotype=meta2.general.animal.genotype;
    if isfield(meta2.imaging.microscope_config,'lightsheet')
        meta.F=meta2.imaging.microscope_config.lightsheet.scanning.z.frequency;
        meta.num_planes = meta2.imaging.microscope_config.lightsheet.scanning.triggering.n_planes;
        meta.piezo_z_amp = meta2.imaging.microscope_config.lightsheet.scanning.z.piezo_max - meta2.imaging.microscope_config.lightsheet.scanning.z.piezo_min;
    else
        meta.F=meta2.imaging.microscope_config.piezo_z.frequency;
        meta.num_planes = meta2.imaging.microscope_config.camera_trigger.n_planes;
        meta.piezo_z_amp = meta2.imaging.microscope_config.piezo_z.amplitude;
    end   
    meta.pulse_times = (0:(1/meta.F)/meta.num_planes:1/meta.F)';
    meta.pulse_times = meta.pulse_times(1:end-1);
    meta.z_step = meta.piezo_z_amp*2/meta.num_planes;
    meta.rez = 0.6;
    if isfield(meta2.stimulus.log{3},'lag')
        meta.lag_condition=meta2.stimulus.log{3}.lag;
    else
        meta.lag_condition=meta2.stimulus.log{4}.lag;
    end
    %% work with stimulus
    filename_stim=cell2mat(dm_dir('*_stimulus_log.hdf5'));
    var_names_stim = h5read(filename_stim,'/data/block0_items');
    for i=1:length(var_names_stim)
        this_name = var_names_stim{i};
        this_name(double(this_name)==0)=[];
        var_names_stim{i} = this_name;
    end
    stim = h5read(filename_stim,'/data/block0_values');
    bad_points=diff(stim(strcmp(var_names_stim,'t'),:))==0;
    stim(:,bad_points)=[];
    t_stim=stim(strcmp(var_names_stim,'t'),:);
    base_grspeed=uint8(-nansum(stim(contains(var_names_stim,'base_vel'),:),1));
    grspeed=single(-nansum(stim(contains(var_names_stim,'vel') & ~contains(var_names_stim,'base_vel'),:),1));
    swim=nansum(stim(contains(var_names_stim,'fish_swimming'),:),1);
    clear stim filename_stim var_names_stim;
    bout_starts_ids=find(diff([0 swim])==1);
    bout_ends_ids=find(diff([0 swim 0])==-1);
    if ~isempty(bout_ends_ids)
        if bout_ends_ids(end)>length(t_stim)
            bout_ends_ids(end)=length(t_stim);
        end
    end
    nob=length(bout_starts_ids);
    
    %% work with behavior
    filename_be=cell2mat(dm_dir('*_behavior_log.hdf5'));
    var_names_be = h5read(filename_be,'/data/block0_items');
    for i=1:length(var_names_be)
        this_name = var_names_be{i};
        this_name(double(this_name)==0)=[];
        var_names_be{i} = this_name;
    end
    be = h5read(filename_be,'/data/block0_values');
    bad_points=diff(be(strcmp(var_names_be,'t'),:))==0;
    be(:,bad_points)=[];
    t_be=be(strcmp(var_names_be,'t'),:);
    tail=be(strcmp(var_names_be,'tail_sum'),:);
    clear be filename_be var_names_be;
    % interpolate tail trace, swim and grating speed
    tail=interp1(t_be,tail,time_be);
    swim=interp1(t_stim,swim,time_be)>=0.5;
    grspeed=interp1(t_stim,grspeed,time_be);
    grmov=interp1(t_stim,double(base_grspeed),time_be)>0.5;
    % scale the tail
    tail=tail-nanmean(tail(~swim));
    tail=tail/nanstd(tail(swim));
    tail(isnan(tail))=0;
    % build vigor trace
    vigor=nan(1,length(tail));
    for i=0.05/dt:length(tail)
        vigor(i)=std(tail(i-(0.05/dt-1):i));
    end
    
    %% detect bouts accurately (and find min and max tail positions)
    bouts.start=nan(1,nob);
    bouts.end=nan(1,nob);
    bouts.bad_bouts=false(1,nob);
    mm_time_array=nan(nob,1000);
    mm_val_array=mm_time_array;
    max_tf_array=false(nob,1000);
    for b=1:nob
        bst0=round((t_stim(bout_starts_ids(b))-0.1)/dt);
        bet0=round(t_stim(bout_ends_ids(b))/dt);
        bouts.start(b)=t_stim(bout_starts_ids(b));
        bouts.end(b)=t_stim(bout_ends_ids(b));
        
        % find individual flicks
        this_mm_time_array=[];
        this_mm_val_array=[];
        this_max_tf_array=[];
        c=0;
        if bst0<2
            bst0=2;
        end
        if bet0>exp_dur/dt-1
            bet0=exp_dur/dt-1;
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
    bouts.main_bouts=t_stim(bout_ends_ids)<t_post_start & t_stim(bout_starts_ids)>=t_main_start;
    bouts.main_bouts=bouts.main_bouts;
    % good bouts: i.e. they are not bad, not spontaneous and happened during main part
    bouts.good_bouts=bouts.main_bouts & ~bouts.bad_bouts & ~bouts.spont_bouts;
    
    %% find bout parameters
    bouts.trial=nan(1,nob);
    bouts.bout_duration=bouts.end-bouts.start;
    bouts.next_interbout_duration=nan(1,nob);
    for g=1:notr
        bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
        for b=bouts_in_this_trial
            bouts.trial(b)=g;
            if b<bouts_in_this_trial(end)
                bouts.next_interbout_duration(b)=bouts.start(b+1)-bouts.end(b);
            end
            bs=round((bouts.start(b)-0.1)/dt);
            be=round(bouts.end(b)/dt);
            if be-bs>=1.1/dt
                be=bs+1.1/dt-1;
            end
            bouts.tail(b,:)=zeros(1,1.1/dt);
            bouts.tail(b,1:be-bs+1)=tail(bs:be);
            bouts.power(b,:)=zeros(1,1.1/dt);
            bouts.power(b,1:be-bs+1)=tail(bs:be).^2;
            bouts.grspeed(b,:)=nan(1,1.1/dt);
            bouts.grspeed(b,1:be-bs+1)=grspeed(bs:be);
        end
    end
    for p=1:length(bout_par_names)
        bouts.(bout_par_names{p})(bouts.bad_bouts | bouts.spont_bouts)=nan;
    end
    for p=1:length(bout_arrays_names)
        bouts.(bout_arrays_names{p})(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    end
    bouts.tail(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
    
    %% find trial averages
    trials=[];
    vigor_for_imagesc=nan(notr,full_trial_dur/dt);
    for g=1:notr
        trial_start=round((gr_starts(g)-rest_dur)/dt)+1;
        trial_end=round((gr_ends(g)+rest_dur)/dt);
        vigor_for_imagesc(g,:)=vigor(trial_start:trial_end);
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
            else
                trials.(this_par).first(g)=nan;
                trials.(this_par).mean(g)=nan;
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
            else
                trials.(this_par).first(g,:)=nan(1,size(bouts.(this_par),2));
                trials.(this_par).mean(g,:)=nan(1,size(bouts.(this_par),2));
            end
        end
    end
    
    %% save the data
    save(filename,'time_be','tail','grspeed','grmov','swim','bouts','trials','meta');
end