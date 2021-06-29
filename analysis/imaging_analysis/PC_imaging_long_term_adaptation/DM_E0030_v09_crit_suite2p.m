% this program compute 4 criteria from scores:
% Criterion 1: difference between trials 11:20 and 21:30
% Criterion 2: difference between trials 21:30 and 61:70
% Criterion 3: difference between trials 61:70 and 71:80
% Criterion 4: difference between trials 71:80 and 111:120
% To compute the difference, we average scores within a block of 10 trials
% and subtract
% To binarize the scores, we do this 1000 times on randomly shuffled
% scores, and check whether the actual difference is within 95 % of the
% null-distribution

%% pathnames
clc; close all; clear all;
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_PC_imaging_results=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_scores=[pathname_PC_imaging_results 'scores\'];
pathname_crit=[pathname_PC_imaging_results 'criteria\'];

%% loop through all fish
cd(pathname_scores);
all_fish=dm_dir('*_f*_scores.mat');
all_fish=strrep(all_fish,'_scores.mat','');
n_fish=length(all_fish);
p_thresh = 0.05;
fish_id = all_fish{1};
load([pathname_scores fish_id '_scores.mat'],'scores_be','scores_im','scores_motor_regr');
n_boots = 10000;
progressbar('fish...','bootstraps...');
for f=1:n_fish
    fish_id = all_fish{f};
    filename = [pathname_crit fish_id '_criteria.mat'];
    load([pathname_scores fish_id '_scores.mat'],'scores_im','scores_be','scores_motor_regr');
    n_rois = size(scores_im,1);
    
    %% compute criteria
    crit_im = nan(n_rois,4);
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
    
    %% compute null-distribution
    H0_im=nan(n_rois,n_boots);
    H0_motor_regr=nan(1,n_boots);
    H0_be=nan(1,n_boots);
    for i=1:n_boots
        rand_ind = [1:10 randperm(110)+10];
        
        scores_im_shuf=scores_im(:,rand_ind);
        temp_block2 = nanmean(scores_im_shuf(:,11:20),2);
        temp_block3 = nanmean(scores_im_shuf(:,21:30),2);
        H0_im(:,i)=temp_block3-temp_block2;
        scores_motor_regr_shuf=scores_motor_regr(rand_ind);
        temp_block2 = nanmean(scores_motor_regr_shuf(11:20));
        temp_block3 = nanmean(scores_motor_regr_shuf(21:30));
        H0_motor_regr(i)=squeeze(temp_block3-temp_block2);
        
        scores_be_shuf = scores_be(rand_ind);
        temp_block2 = nanmean(scores_be_shuf(11:20));
        temp_block3 = nanmean(scores_be_shuf(21:30));
        H0_be(i)=temp_block3-temp_block2;
        progressbar([],i/n_boots);
    end
    
    %% find percentiles
    prcnt_min_im = prctile(H0_im,p_thresh*100/2,2);
    prcnt_max_im = prctile(H0_im,100-p_thresh*100/2,2);
    prcnt_min_motor_regr = prctile(H0_motor_regr,p_thresh*100/2);
    prcnt_max_motor_regr = prctile(H0_motor_regr,100-p_thresh*100/2);
    prcnt_min_be = prctile(H0_be,p_thresh*100/2);
    prcnt_max_be = prctile(H0_be,100-p_thresh*100/2);
    
    %% find significant criteria
    crit_im_signif = zeros(n_rois,4,'single');
    for i=1:4
        crit_im_signif(crit_im(:,i)<prcnt_min_im,i)=-1;
        crit_im_signif(crit_im(:,i)>prcnt_max_im,i)=1;
    end
    crit_motor_regr_signif = zeros(1,4,'single');
    crit_motor_regr_signif(crit_motor_regr<prcnt_min_motor_regr)=-1;
    crit_motor_regr_signif(crit_motor_regr>prcnt_max_motor_regr)=1;
    crit_be_signif = zeros(1,4);
    crit_be_signif(crit_be<prcnt_min_be)=-1;
    crit_be_signif(crit_be>prcnt_max_be)=1;
    
    %% save
    save(filename,'crit_be','crit_be_signif','crit_im','crit_im_signif','crit_motor_regr','crit_motor_regr_signif','H0_im','H0_motor_regr','H0_be'); 
    progressbar(f/n_fish,[]);
end