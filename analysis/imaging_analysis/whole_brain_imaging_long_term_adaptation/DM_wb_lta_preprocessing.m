% this progam performs pre-processing of the traces and ROI
% of the long-term adaptation experiment during whole-brain lightsheet imaging

%% pathnames
clc; close all; clear all;
% path to raw data
pathname_MS_data = 'C:\Markov_et_al_2021_Nat_Commun_data_code\data\';
pathname_wb_imaging_results=[pathname_MS_data 'whole_brain_imaging_long_term_adaptation\'];
pathname_raw_ROIs = [pathname_wb_imaging_results 'raw_ROIs\'];
pathname_raw_traces = [pathname_wb_imaging_results 'raw_traces\'];
pathname_affine_x=[pathname_wb_imaging_results 'anatomies\morphed\'];
pathname_temp_lists='C:\users\dmark\';
pathname_ROIs = [pathname_wb_imaging_results 'ROIs\'];
pathname_traces = [pathname_wb_imaging_results 'traces\'];
pathname_behavior=[pathname_wb_imaging_results 'behavior\'];

%% size and rezolution of the reference brain
[sz_ref,rez_ref]=read_nrrd_metadata([pathname_MS_data 'reference_brain_stacks\PortuguesLab_wholebrain_ref.nrrd']);
% whole brian mask
ref_mask=nrrdread([pathname_MS_data 'reference_brain_stacks\PortuguesLab_wholebrain_ref_mask.nrrd']);
ref_mask=ref_mask>0;

%% cutoff frequency for high-pass filtering
fc_baseline=1/300; % Hz, 5 min period

%% for building the motor regressor
s_half_decay = genotype2tau('GCaMP6s');

%% loop through all fish
all_fish=dm_dir([pathname_raw_traces '*_f*']);
n_fish=length(all_fish);
progressbar('Fish progress...','Morphing ROIs...','Fixing ROI IDs, processing traces...');
for f=18%1:n_fish
    fish_id=all_fish{f};
    
    %% load traces
    traces = h5read([pathname_raw_traces fish_id '\traces.h5'],'/traces')';
    traces=double(traces);
    traces(all(isnan(traces),2),:)=[];
    [n_ROIs, n_frames] = size(traces);
    
    %% load ROIs
    ROIs_local = tiff2mat(pathname_raw_ROIs, ['rois_' fish_id '.tif']);
    ROIs_local = ROIs_local - min(ROIs_local(:))+1;
    ROIs_local(ROIs_local==max(ROIs_local(:))) = 0;
    ROIs_local=flip(ROIs_local,3);
    ROIs_local = rot90(ROIs_local);
    sz_ROIs_local=size(ROIs_local);
    
    %% affine x folder for morphing ROIs
    affine_x=[pathname_affine_x fish_id '_anatomy_affine.xform'];
    
    %% behavior
    load([pathname_behavior fish_id '_behavior.mat'],'meta','swim','time_be');
    pulse_times=flip(meta.pulse_times);
    
    %% high-pass filter for removing slow drift in the traces
    [bb_baseline,aa_baseline]=butter(3,2*fc_baseline/meta.F,'high');
    traces(traces==0)=nan;
    
    %% process this fish
    % get IDis of remaining ROIs
    all_ROIs=unique(ROIs_local)';
    all_ROIs(all_ROIs==0)=[];  
    % create output variables
    ROIs=zeros(sz_ref,'uint32');
    time_offsets=nan(n_ROIs,1);
    % start the loop that morphs the ROIs to the refernce brain
    c=0;
    for i=all_ROIs
        c=c+1;
        coord0=find(ROIs_local==i);
        [row0,column0,plane0]=ind2sub(sz_ROIs_local, coord0);
        ori_plane=plane0(1);
        row0=row0*meta.rez;
        column0=column0*meta.rez;
        plane0=plane0*meta.z_step;
        ori_plane_um=ori_plane*meta.z_step;
        extended_planes=ori_plane_um-1:-1:ori_plane_um-meta.z_step+1;
        num_ext_planes=length(extended_planes);
        for ii=extended_planes
            plane0=[plane0; ones(length(row0),1)*ii];
        end
        row0=repmat(row0,num_ext_planes+1,1);
        column0=repmat(column0,num_ext_planes+1,1);
        
        coord1=DM_morph_non_ref_coord([row0,column0,plane0],affine_x,pathname_temp_lists);
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
                % find pixels that are intesecting and remove 1/2
                px = find(this_ROI & ROIs==ii);
                this_ROI(px(randperm(length(px),round(length(px)/2))))=false;
            end
            ROIs(this_ROI)=i;
            % write the temporal offset of this roi
            time_offsets(i)=pulse_times(ori_plane);
        end
        progressbar([],c/n_ROIs,[])
    end
    % remove out of brain ROIs
    ROIs(~ref_mask)=0;
    % find all remaining ROIs
    all_ROIs=unique(ROIs)';
    all_ROIs(all_ROIs==0)=[];
    n_ROIs=length(all_ROIs);
    % fix ROI IDs and process the traces
    ROIs2=zeros(sz_ref,'uint32');
    traces2=nan(n_ROIs,n_frames,'single');
    time_offsets2=nan(n_ROIs,1);
    ROI_coord=cell(n_ROIs,1);
    for i=1:n_ROIs
        old_id=all_ROIs(i);
        coord = ROIs==old_id;
        ROIs2(coord)=i;
        time_offsets2(i)=time_offsets(old_id);
        traces2(i,:) = process_trace_here(traces(old_id,:),bb_baseline,aa_baseline);   
        ROI_coord{i}=int32(find(coord));
        progressbar([],[],i/n_ROIs);
    end
    ROIs=ROIs2;
    traces=traces2;
    time_offsets=time_offsets2;
    
    %% build motor regressor
    time_im=1/meta.F:1/meta.F:n_frames/meta.F;
    my_kernel=2.^(-time_be/s_half_decay);
    trace_motor_regr=trace2regressor(swim,my_kernel,time_be,time_im);
    trace_motor_regr = process_trace_here(trace_motor_regr,bb_baseline,aa_baseline);  
 
    %% save
    save([pathname_ROIs fish_id '_ROIs.mat'],'ROIs','ROI_coord');
    save([pathname_traces fish_id '_traces.mat'],'traces','time_offsets','trace_motor_regr');
    progressbar(f/n_fish,[],[])
end

function [trace] = process_trace_here(trace,bb_baseline,aa_baseline)
my_nans=isnan(trace);
trace(my_nans)=nanmean(trace);
trace=filtfilt(bb_baseline,aa_baseline,trace);
trace(my_nans)=nan;
trace=trace-nanmean(trace);
trace=trace/nanstd(trace);
trace=single(trace);
end