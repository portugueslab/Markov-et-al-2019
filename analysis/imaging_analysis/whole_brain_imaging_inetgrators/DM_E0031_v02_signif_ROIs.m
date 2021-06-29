%% pathnames
clc; close all; clear all;
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_ROIs=[pathname_imaging_results 'ROIs\'];
pathname_clustering=[pathname_imaging_results 'clustering\'];
pathname_time_constants=[pathname_imaging_results 'time_constants\'];
pathname_reference_brains = [pathname_MS_data 'reference_brain_stacks\'];

%% size of the reference brain
sz=read_nrrd_metadata([pathname_reference_brains 'PortuguesLab_wholebrain_ref.nrrd']);

%% load the data (coordinates and time constants)
all_fish=dm_dir([pathname_ROIs '*_f*_ROIs.mat']);
all_fish=strrep(all_fish,'_ROIs.mat','');
n_fish=length(all_fish);
coord={};
clust=uint8([]);
for f=1:n_fish
    fish_id=all_fish{f};
    load([pathname_ROIs fish_id '_ROIs.mat'],'ROI_coord');
    coord=[coord; ROI_coord];
    load([pathname_clustering fish_id '_clustering.mat'],'sensmot_clust');
    load([pathname_time_constants fish_id '_time_constants.mat'],'time_constants');
    sensmot_clust(time_constants>1.5)=3;
    clust=[clust; sensmot_clust];
end
n_ROIs = length(clust);
clear time_constants ROI_coord f all_fish n_fish fish_id sensmot_clust;

%% build stacks where value of each voxel equals number of fish with activity in that voxel
% for all active ROIs and for individual functional labels
stack_num_fish=zeros(sz,'uint8');
stack_num_sensors=stack_num_fish;
stack_num_motor=stack_num_fish;
stack_num_integrators=stack_num_fish;
for i=1:n_ROIs
    this_coord=coord{i};
    switch clust(i)
        case 1
            stack_num_fish(this_coord)=stack_num_fish(this_coord)+1;
            stack_num_sensors(this_coord)=stack_num_sensors(this_coord)+1;
        case 2
            stack_num_fish(this_coord)=stack_num_fish(this_coord)+1;
            stack_num_motor(this_coord)=stack_num_motor(this_coord)+1;
        case 3
            stack_num_fish(this_coord)=stack_num_fish(this_coord)+1;
            stack_num_integrators(this_coord)=stack_num_integrators(this_coord)+1;
    end
end
clear i this_coord;

%% Formulate the H0
% All active ROIs have equal probability to be asiged with any of the three
% funtional labels (p=1/3).
% Note: active ROI means that it has any of the functional labels.

% We test this H0 for each ROI. To do so, we first define a term "locus"
% as all voxels that belong to that ROI. We then find number of
% fish with active ROIs in that locus (N) and number of fish with ROIs that
% share the functional label with the original ROI in that locus (M). We then
% estimate the probability P of observing M or more ROIs out of total N active ROIs given H0.
% if P < 0.05, we reject the H0 and conclude that ROIs in that locus are
% more likely to belong to the functional group of the original ROI than to
% any other functional group. 
p=1/3;
signif_ROIs=false(n_ROIs,1);
p_thresh = 0.05;
for i=1:n_ROIs
    if clust(i)~=0 % if this ROI is active
        % find its "locus"
        this_coord=coord{i};
        % find number of fish with active ROIs in that locus
        N = max(stack_num_fish(this_coord));
        % find number of fish with ROIs with the same functional label in that locus
        switch clust(i)
            case 1
                M=max(stack_num_sensors(this_coord));
            case 2
                M=max(stack_num_motor(this_coord));
            case 3
                M=max(stack_num_integrators(this_coord));
        end
        % compute the probability P of this or more extreme observations given the H0
        % (by "more extreme" I mean observing M or more ROIs)
        P=0;
        for ii=M:N
            P=P+binopdf(ii,N,1/3);
        end
        % make a decision about rejecting the H0 given 95% significance level
        if P<p_thresh
           signif_ROIs(i)=true;
        end
    end
end
save([pathname_imaging_results 'significant_ROIs.mat'],'signif_ROIs');