%% initial stuff
clc; close all; clear all;
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\imaging_analysis\PC_imaging_long_term_adaptation');
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\MATLAB_functions');

%% morphing anatomies
% Manually save _anatomy.tif as .nrrd using FIJI nrrd writer
% Danger! Do not run the next program as a whole
% Morphing of these data was tricky, so
% it should be run section by section with manual tweaking
DM_E0030_v09_morph_anatomies;

%% pre-processing
% analysis of the behavior, same as for behavioral experiments
DM_E0030_v09_behavior;
% extract traces and ROIs, filters traces, build motor regressor
DM_E0030_v09_preprocessing_suite2p;
% morph ROIs to the reference brain
DM_E0030_v09_morph_ROIs;

% these four programs save pre-processed data to:
% ...\Markov_et_al_2021_Nat_Commun_data&code\data\PC_imaging_long_term_adaptation\...
% behavior
% anatomies 
% traces
% ROIs

%% subsequent analysis
% compute grating- and bout-triggered averages
DM_E0030_v09_trigaver_suite2p;
% compte scores
DM_E0030_v09_scores_suite2p;
% compute criteria
DM_E0030_v09_crit_suite2p;

% these three programs save final data to:
% ...\Markov_et_al_2021_Nat_Commun_data&code\data\PC_imaging_long_term_adaptation\...
% triggered_traces
% scores
% criteria