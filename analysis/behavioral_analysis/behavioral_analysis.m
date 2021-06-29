% this program analyzes the raw data (detects individual bouts and computes stuff used in the figures)
clc; close all; clear all;
addpath('...\MATLAB_functions');

%% pathnames
pathname_MS_data = '...\data\';
pathname_behavior=[pathname_MS_data 'behavior\'];
exp_names={'acute_reaction_experiment\',...
    'long_term_adaptation_experiment\normal_reafference_control\',...
    'long_term_adaptation_experiment\lag_trained\'};
exp_groups={'WT_TL_group','treatment_control_group','PC_ablated_group'};

%% timing
dt=0.005; % [s]
num_trials=240;
trial_dur=30; % [s]
time_be=dt:dt:num_trials*trial_dur;
num_frames=length(time_be);
vigor_window_frames=round(0.05/dt);
gr_starts=7.5:30:7.5+239*30;
gr_ends=gr_starts+15;
t_main_start=20*trial_dur;
t_post_start=230*trial_dur;

%% variables with parameter and reafference condition names
bout_par_names={'bout_duration','next_interbout_duration','power'};
cond.names={'lag';'shunted_lag';'gain';'gain_drop'};
cond.xlabel=[{{'0','75','150','225','300','Inf'}};...
    {{'0','75','150','225','300','Inf'}};...
    {{'0', '0.33','0.66','1','1.33','1.66','2'}};...
    {{'1111','0111','0011','0001','0000','1110','1100','1000'}}];
cond.vals=[{[0 0.075 0.15 0.225 0.3 0; 1 1 1 1 1 0; 0 0 0 0 0 0; 0 0 0 0 0 0]};...
    {[0 0 0 0 0 0; 1 1 1 1 1 0; 0 0 0 0 0 0; 0 0.075 0.15 0.225 0.3 0]};...
    {[0 0 0 0 0 0 0; 0 0.33 0.66 1 1.33 1.66 2; 0 0 0 0 0 0 0; 0 0 0 0 0 0 0]};...
    {[0 0 0 0 0 0 0 0; 1 1 1 1 1 1 1 1; 0 0 0 0 0 0.225 0.15 0.075; 0 0.075 0.15 0.225 0.3 0.3 0.3 0.3]}];

%% main loop
progressbar('Experiment...','Genotype...','Fish...');
for exp=1:length(exp_names)
    for gen=1:length(exp_groups)
        behavior=struct;
        pathname=[pathname_behavior exp_names{exp} exp_groups{gen}];
        cd(pathname);
        all_fish=dm_dir('*_f*_data.mat');
        for f=1:length(all_fish)
            
            %% load data from this fish
            load(all_fish{f},'tail','time_tail','time_stim','swim','metadata');
            if contains(exp_names{exp},'acute_adaptation')
                load(all_fish{f},'lag','gain','gd_starts','gd_ends');
            end
            
            %% find bout starts and ends (detected online during the experiment)
            swim(1)=false; swim(end)=false;
            bout_starts_ids=find(diff([false swim])==1);
            bout_ends_ids=find(diff([false swim false])==-1);
            num_bouts=length(bout_starts_ids);
            bouts=[];
            
            %% find reafference of all these bouts
            % (if this is acute adaptation experiment)
            if contains(exp_names{exp},'acute_adaptation')
                bouts.gain=gain(bout_starts_ids)';
                bouts.lag=lag(bout_starts_ids)';
                bouts.gain_drop_start=gd_starts(bout_starts_ids)';
                bouts.gain_drop_start(isnan(bouts.gain_drop_start))=0;
                bouts.gain_drop_end=gd_ends(bout_starts_ids)';
                bouts.gain_drop_end(isnan(bouts.gain_drop_end))=0;
            end
            
            %% work with behavioral traces
            % interpolate tail and swim to time_be
            tail=interp1(time_tail,tail,time_be);
            swim=interp1(time_stim,double(swim),time_be)>0;
            % scale the tail trace
            tail=tail-nanmean(tail(~swim));
            tail=tail/nanstd(tail(swim));
            tail(isnan(tail))=0;
            % build vigor trace
            vigor=nan(1,num_frames);
            for i=vigor_window_frames:num_frames
                vigor(i)=std(tail(i-(vigor_window_frames-1):i));
            end
             
             %% detect bouts accurately by identifying individual tail flicks
             bouts.start=nan(num_bouts,1,'single');
             bouts.end=nan(num_bouts,1,'single');
             bouts.bad_bouts=false(num_bouts,1);
             mm_time_array=nan(num_bouts,1000);
             mm_val_array=mm_time_array;
             max_tf_array=false(num_bouts,1000);
             for b=1:num_bouts
                 bst0=round((time_stim(bout_starts_ids(b))-0.1)/dt);
                 bet0=round(time_stim(bout_ends_ids(b))/dt);
                 bouts.start(b)=time_stim(bout_starts_ids(b));
                 bouts.end(b)=time_stim(bout_ends_ids(b));
                 
                 % find individual flicks
                 this_mm_time_array=[];
                 this_mm_val_array=[];
                 this_max_tf_array=[];
                 c=0;
                 if bst0<2
                     bst0=2;
                 end
                 if bet0>240*30/dt-1
                     bet0=240*30/dt-1;
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
                 % remove all "bad" flicks
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
                 [false; bouts.start(2:end)-bouts.end(1:end-1)<0.1] | [bouts.start(2:end)-bouts.end(1:end-1)<0.1; false] |... % bouts with interbouts shorter than 100 ms (tipically, these are bouts which were detected as 2 bouts)
                 all(isnan(mm_time_array),2) |...  % wierd bouts with no flicks
                 max(diff(mm_time_array,1,2),[],2)>0.1; % bouts with max delta flick time > 100 ms (this happens during some wierd ugly bouts)
             
             % find spontaneous bouts (i.e. bouts which started after trial start and finished before trial end)
             bouts.spont_bouts=true(num_bouts,1);
             for g=1:num_trials
                 bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
                 bouts.spont_bouts(bouts_in_this_trial)=false;
             end
             % find bouts during the main part of the experiment (adaptation phase)
             bouts.main_bouts=time_stim(bout_ends_ids)<t_post_start & time_stim(bout_starts_ids)>=t_main_start;
             bouts.main_bouts=bouts.main_bouts';
             % good bouts: i.e. they are not bad, not spontaneous and happened during main part
             bouts.good_bouts=bouts.main_bouts & ~bouts.bad_bouts & ~bouts.spont_bouts;
            
             %% find bout parameters
             bouts.trial=zeros(num_bouts,1,'uint8');
             bouts.bout_duration=bouts.end-bouts.start;
             bouts.next_interbout_duration=nan(num_bouts,1,'single');
             for g=1:num_trials
                 bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
                 for b=bouts_in_this_trial'
                     bouts.trial(b)=g;
                     if b<bouts_in_this_trial(end)
                         bouts.next_interbout_duration(b)=bouts.start(b+1)-bouts.end(b);
                     end
                     bs=round((bouts.start(b)-0.1)/dt);
                     be=round(bouts.end(b)/dt);
                     if be-bs>=1.1/dt
                         be=bs+1.1/dt-1;
                     end
                     bouts.power(b,:)=zeros(1,1.1/dt,'single');
                     temp_tail=tail(bs:be);
                     temp_tail=temp_tail-nanmedian(temp_tail(1:0.1/dt));
                     bouts.power(b,1:be-bs+1)=temp_tail.^2;
                 end
             end
             for p=1:length(bout_par_names)
                 bouts.(bout_par_names{p})(bouts.bad_bouts | bouts.spont_bouts,:)=nan;
             end
            
             %% find trial averages
             trials=[];
             for g=1:num_trials
                 trial_start=round((gr_starts(g)-7.5)/dt)+1;
                 trial_end=round((gr_ends(g)+7.5)/dt);
                 bouts_in_this_trial=find(bouts.start>gr_starts(g) & bouts.end<gr_ends(g));
                 trials.bouts_in_this_trial(g,1)={bouts_in_this_trial};
                 trials.number_of_bouts(g,1)=length(bouts_in_this_trial);
                 if ~isempty(bouts_in_this_trial)
                     trials.latency(g,1)=bouts.start(bouts_in_this_trial(1))-gr_starts(g);
                 else
                     trials.latency(g,1)=nan;
                 end
                 for p=1:length(bout_par_names)
                     this_par=bout_par_names{p};
                     these_bouts=bouts.(this_par)(bouts_in_this_trial,:);
                     these_bouts=these_bouts(~all(isnan(these_bouts),2),:);
                     trials.(this_par).n(g,1)=size(these_bouts,1);
                     if trials.(this_par).n(g)>0
                         trials.(this_par).first(g,:)=these_bouts(1,:);
                         trials.(this_par).mean(g,:)=nanmean(these_bouts,1);
                     else
                         trials.(this_par).first(g,:)=nan(1,size(bouts.(this_par),2));
                         trials.(this_par).mean(g,:)=nan(1,size(bouts.(this_par),2));
                     end
                 end
             end
             
             %% find condition averages (for acute reaction experiment)
             if contains(exp_names{exp},'acute_adaptation')
                 for c=1:length(cond.names)
                     this_cond_name=cond.names{c};
                     this_cond_vals=single(cond.vals{c});
                     bouts.condition_mean.(this_cond_name).conditions=cond.xlabel{c};
                     for i=1:length(bouts.condition_mean.(this_cond_name).conditions)
                         tf=bouts.good_bouts & bouts.lag==this_cond_vals(1,i) & bouts.gain==this_cond_vals(2,i) & bouts.gain_drop_start==this_cond_vals(3,i) & bouts.gain_drop_end==this_cond_vals(4,i);
                         for p=1:length(bout_par_names)
                             this_par=bout_par_names{p};
                             these_bouts=bouts.(this_par)(tf,:);
                             these_bouts=these_bouts(~all(isnan(these_bouts),2),:);
                             bouts.condition_mean.(this_cond_name).(this_par).n(i)=size(these_bouts,1);
                             if bouts.condition_mean.(this_cond_name).(this_par).n(i)>0
                                 bouts.condition_mean.(this_cond_name).(this_par).mean(i,:)=nanmean(these_bouts,1);
                             else
                                 bouts.condition_mean.(this_cond_name).(this_par).mean(i,:)=nan(1,size(bouts.(this_par),2));
                             end
                         end
                     end
                 end
             end
             
             %% save pooled data into the final structure
             behavior(f).fish_id=strrep(all_fish{f},'_data.mat','');
             behavior(f).metadata=metadata;
             behavior(f).bouts=bouts;
             behavior(f).trials=trials;
             progressbar([],[],f/length(all_fish));            
        end
        save('pooled_data.mat','behavior');
        progressbar([],gen/length(exp_groups),[]);  
    end
    progressbar(exp/length(exp_names),[],[]); 
end
    