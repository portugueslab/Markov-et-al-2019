%% initial stuff
clc; close all; clear all;
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\imaging_analysis\whole_brain_imaging_inetgrators');
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\MATLAB_functions');

%% morphing anatomies
% Manually save _anatomy.tif as .nrrd using FIJI nrrd writer
% Danger! Do not run the next program as a whole
% Morphing of these data was tricky, so
% it should be run section by section with manual tweaking
% DM_E0031_v02_morph_anatomies;

%% pre-processing
% reorganize the raw data and save as mat files
% DM_E0031_v02_compress_csv_files;
% analysis of the behavior, same as for behavioral experiments
% DM_E0031_v02_behavior;
% pre-process traces, morph ROIs
% DM_E0031_v02_preprocessing;

% these four programs save pre-processed data to:
% ...\data\whole_brain_imaging_inetgrators\...
% behavior
% anatomies 
% traces
% ROIs

%% subsequent analysis
% compute triggered responses, scores, and cluster ROIs into sensory and motor
DM_E0031_v02_trigaver_scores_clustering;
% fit time constnts to sensory ROIs and cluster them into sensors and integrators
DM_E0031_v02_fit_time_constants;
% find location of ROIs that are consistent across fish
DM_E0031_v02_signif_ROIs;