% this will morph ROIs, compute their coordinates and pre-process traces

%% pathnames
clc; close all; clear all;
% path to raw data
pathname_raw_data='J:\_Shared\experiments\E0031_acute_adaptation\v02_lightsheet\data_imaging\';
path_to_affine_x='C:\Users\dmarkov\Desktop\Data_Daniil\Internal_models_figures\morphing\flipped_anatomies_morphed\';
path_to_lists='C:\Users\dmarkov\Desktop\Data_Daniil\Internal_models_figures\morphing\temp_pathname\';
% path to pre-processed data uploaded with the manuscript
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_imaging_results=[pathname_MS_data 'whole_brain_imaging_inetgrators\'];
pathname_ROIs=[pathname_imaging_results 'ROIs\'];
pathname_traces=[pathname_imaging_results 'traces\'];
pathname_anatomies = [pathname_imaging_results 'anatomies\'];
pathname_behavior = [pathname_imaging_results 'behavior\'];

%% size and rezolution of the reference brain
[sz_ref,rez_ref]=read_nrrd_metadata('C:\Users\dmark\OneDrive\Рабочий стол\internal_models\Markov_et_al_2021_Nat_Commun_data&code\data\reference_brain_stacks\PortuguesLab_wholebrain_ref.nrrd');
% whole brian mask
ref_mask=nrrdread('C:\Users\dmark\OneDrive\Рабочий стол\internal_models\Markov_et_al_2021_Nat_Commun_data&code\data\reference_brain_stacks\PortuguesLab_wholebrain_ref_mask.nrrd');
ref_mask=ref_mask>0;

%% cutoff frequency for high-pass filtering
fc_baseline=1/300; % Hz, 5 min period

%% loop through all fish
all_fish=dm_dir([pathname_anatomies '*_f*_anatomy.nrrd']);
all_fish=strrep(all_fish,'_anatomy.nrrd','');
n_fish=length(all_fish);
progressbar('Fish progress...','Morphing ROIs...','Fixing ROI IDs, processing traces...','Computing ROI coordinates...')
for f=1:n_fish
    this_name=all_fish{f};
    
    %% load traces, ROIs, as well as  metadata
    % path to raw data of this fish
    path_to_fish=[pathname_raw_data this_name '\exp0031\'];
    affine_x=[path_to_affine_x this_name '_anatomy.xform'];
    % load ROIs and traces, metadata, prepare time_offsets
    ROIs_local = tiff2mat(path_to_fish, 'rois.tif')+1;
    sz_ROIs_local=size(ROIs_local);
    traces=load([path_to_fish 'matlab_data.mat'],'traces');
    traces=traces.traces';
    meta=load([pathname_behavior this_name '_behavior.mat'],'meta');
    meta=meta.meta;
    pulse_times=flip(meta.pulse_times);
    
    %% filters for computing baseline and denoising the trace
    [bb_baseline,aa_baseline]=butter(3,2*fc_baseline/meta.fs,'low');
    fc=1/genotype2tau(meta.genotype);
    [bb_6s,aa_6s]=butter(3,2*fc/meta.fs,'low');
    
    %% process this fish
    % remove ROIs from the first 4 planes
    ROIs_local(:,:,1:4)=0;
    % remove ROIs with artifact traces
    max_diff=max(abs(diff(traces,1,2)),[],2);
    artifact_ROIs=find(max_diff>100);
    ROIs_local(ismember(ROIs_local,artifact_ROIs))=0;
    % get IDis of remaining ROIs
    all_ROIs=unique(ROIs_local)';
    all_ROIs(all_ROIs==0)=[];
    n_ROIs=length(all_ROIs);
    % flip ROIs
    ROIs_local=flip(ROIs_local,3);
    % create output variables
    ROIs=zeros(sz_ref,'uint32');
    time_offsets=nan(n_ROIs,1);
    % start the loop that morphs the ROIs to the refernce brain
    c=0;
    for i=all_ROIs
        c=c+1;
        coord0=find(ROIs_local==i);
        [row0,column0,plane0]=ind2sub(sz_ROIs_local, coord0);
        ori_plane=unique(plane0);
        row0=row0*meta.rez(1);
        column0=column0*meta.rez(2);
        plane0=plane0*meta.rez(3);
        ori_plane_um=ori_plane*meta.rez(3);
        extended_planes=ori_plane_um-1:-1:ori_plane_um-meta.rez(3)+1;
        num_ext_planes=length(extended_planes);
        for ii=extended_planes
            plane0=[plane0; ones(length(row0),1)*ii];
        end
        row0=repmat(row0,num_ext_planes+1,1);
        column0=repmat(column0,num_ext_planes+1,1);
        
        coord1=DM_morph_non_ref_coord([row0,column0,plane0],affine_x,path_to_lists);
        coord1=coord1./rez_ref;
        coord1=unique([floor(coord1);ceil(coord1)],'rows');
        bad_coord=find(any(coord1<=0,2) | coord1(:,1)>sz_ref(1) | coord1(:,2)>sz_ref(2) | coord1(:,3)>sz_ref(3));
        coord1(bad_coord,:)=[];
        if ~isempty(coord1)
            % write it's mid plane
            % add this ROI coordinates
            coord1=sub2ind(sz_ref, coord1(:,1), coord1(:,2), coord1(:,3));
            this_ROI = false(sz_ref);
            this_ROI(coord1)=true;
            % find if there are rois already
            ROIs_there = unique(ROIs (this_ROI))';
            ROIs_there(ROIs_there==0)=[];
            for ii=ROIs_there
                % find pixels that are intesepting and remove 1/2
                px = find(this_ROI & ROIs==ii);
                this_ROI(px(randperm(length(px),round(length(px)/2))))=false;
            end
            ROIs(this_ROI)=i;
            % write the temporal offset of this roi
            time_offsets(i)=pulse_times(ori_plane);
        end
        progressbar([],c/n_ROIs,[],[])
    end
    % remove out of brain ROIs
    ROIs(~ref_mask)=0;
    % find all remaining ROIs
    all_ROIs=unique(ROIs)';
    all_ROIs(all_ROIs==0)=[];
    n_ROIs=length(all_ROIs);
    % fix ROI IDs and process the traces
    ROIs2=zeros(sz_ref,'uint32');
    traces=double(traces);
    traces(traces==0)=nan;
    traces2=nan(n_ROIs,size(traces,2),'single');
    time_offsets2=nan(n_ROIs,1);
    for i=1:n_ROIs
        old_id=all_ROIs(i);
        ROIs2(ROIs==old_id)=i;
        time_offsets2(i)=time_offsets(old_id);
        traces2(i,:) = process_trace_here(traces(old_id,:),bb_baseline,aa_baseline,bb_6s,aa_6s);
        progressbar([],[],i/n_ROIs,[]);
    end
    ROIs=ROIs2;
    traces=traces2;
    time_offsets=time_offsets2;
    % extract linear coordinates of morphed ROIs
    ROI_coord=cell(n_ROIs,1);
    for i=1:n_ROIs
        coord=uint32(find(ROIs==i));
        ROI_coord{i}=coord;
        progressbar([],[],[],i/n_ROIs);
    end
    % save the data
    filename=[pathname_ROIs this_name '_ROIs.mat'];
    save(filename,'ROI_coord');
    filename=[pathname_traces this_name '_traces.mat'];
    save(filename,'traces','time_offsets');
    progressbar(f/length(all_fish),[],[],[])
end

function [trace] = process_trace_here(trace,bb_baseline,aa_baseline,bb_6s,aa_6s)
my_nans=isnan(trace);
trace(my_nans)=nanmean(trace);
trace=filtfilt(bb_6s,aa_6s,trace);
trace_baseline=filtfilt(bb_baseline,aa_baseline,trace);
trace=trace-trace_baseline;
trace(my_nans)=nan;
trace=trace-nanmean(trace);
trace=trace/nanstd(trace);
trace=single(trace);
end