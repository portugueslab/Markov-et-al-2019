% this will compress new stytra acute fish

%% initiate
clc; close all; clear all;

path_to_data='J:\_Shared\experiments\E0031_acute_adaptation\v02_lightsheet\data_behavior\';
cd(path_to_data);
all_fish=dir('*_f*');
%% loop through fish
progressbar('Fish');
for f=1:length(all_fish)
    %% check if this fish is not already compressed
    this_name=all_fish(f).name;
    path_to_fish=[path_to_data this_name];
    cd(path_to_fish);
    filename_data=[this_name '_data.mat'];
    %% compress
    if exist(filename_data,'file')~=2 % only if it was not compressed before
        %% stimulus
        filename_stim=dir('*_stimulus_log.csv');
        A=readtable(filename_stim.name);
        stim=[];
        stim.time=A{:,'t'};
        stim.grspeed_prepost=-A{:,'general_cl1D_vel'};
        stim.base_grspeed_prepost=-A{:,'general_cl1D_base_vel'};
        stim.swim_prepost=A{:,'general_cl1D_fish_swimming'};
        stim.grspeed_main=-A{:,'acute_cl1D_vel'};
        stim.base_grspeed_main=-A{:,'acute_cl1D_base_vel'};
        stim.swim_main=A{:,'acute_cl1D_fish_swimming'};
        stim.gain_main=A{:,'acute_cl1D_gain'};
        stim.lag_main=A{:,'acute_cl1D_lag'};
        stim.gd_starts_main=A{:,'acute_cl1D_gain_drop_start'};
        stim.gd_ends_main=A{:,'acute_cl1D_gain_drop_end'};
        stim.shunted_lag_main=A{:,'acute_cl1D_shunted'};
         
        %% estimator
        filename_estim=dir('*_estimator_log.csv');
        A=readtable(filename_estim.name);
        estim=[];
        estim.time=A{:,'t'};
        estim.vigor=A{:,'vigor'};

        %% behavior
        filename_be=dir('*_behavior_log.csv');
        A=readtable(filename_be.name);
        be=[];
        be.time=A{:,'t'};
        be.tail=A{:,'tail_sum'};
        var_names=A.Properties.VariableNames;
        n_tailpoints=0;
        for i=1:length(var_names)
            if ~isempty(strfind(var_names{i},'theta_'))
                n_tailpoints=n_tailpoints+1;
            end
        end
        for i=1:n_tailpoints
            be.individual_tailpoints(:,i)=A{:,['theta_' num2str(i-1,'%02.f')]};
        end
        
        %% metadata
        filename_meta=dir('*_metadata.json');
        A=read_json(filename_meta.name);
        meta=A;
        
        %% save
        save(filename_data,'stim','estim','be','meta','-v7.3');
        
    end
    progressbar(f/length(all_fish));
end