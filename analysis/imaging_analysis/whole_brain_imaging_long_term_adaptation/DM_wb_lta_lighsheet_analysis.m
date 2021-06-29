%% initial stuff
clc; close all; clear all;
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\imaging_analysis\whole_brain_imaging_long_term_adaptation');
addpath('C:\Markov_et_al_2021_Nat_Commun_data_code\analysis_code\MATLAB_functions');

%% morphing anatomies
% same is in the whole-brain exp

%% pre-processing
% analysis of the behavior, same as for behavioral experiments
DM_wb_lta_behavior;
% pre-process traces, morph ROIs
DM_wb_lta_preprocessing;

% these two programs save pre-processed data to:
% ...\data\whole_brain_imaging_long_term_adaptation\...
% behavior
% anatomies 
% traces
% ROIs

%% subsequent analysis
DM_wb_lta_main_analysis;

% these three programs save final data to:
% ...\data\whole_brain_imaging_long_term_adaptation\processed_data\