%% pathnames
clc; close all; clear all;

% path to raw anatomies
pathname_data='J:\_Shared\experiments\E0030_long_term_adaptation\v09_lightsheet\data_imaging_suite2p\';
pathname_raw_anat=[pathname_data 'anatomies\'];
% path to morphed ROIs
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_PC_imaging_results=[pathname_MS_data 'PC_imaging_long_term_adaptation\'];
pathname_morphed_anat=[pathname_PC_imaging_results 'anatomies\'];
pathname_ROIs = [pathname_PC_imaging_results 'ROIs\'];
path_to_lists=[pathname_PC_imaging_results 'temp_pathname\'];
pathname_ref = [pathname_MS_data 'reference_brain_stacks\'];
filename_ref = 'PortuguesLab_PC_ref.nrrd';
str_ref=[pathname_ref filename_ref];
[sz_ref, rez_ref] = read_nrrd_metadata(str_ref);

%% loop through fish
cd(pathname_ROIs);
all_fish=dm_dir('*_f*_ROIs.mat');
all_fish=strrep(all_fish,'_ROIs.mat','');
n_fish=length(all_fish);
progressbar('Fish progress...','Morphing ROIs...','Computing ROI coordinates...');
for f=1:n_fish
    fish_id=all_fish{f};
    load([pathname_ROIs fish_id '_ROIs.mat'],'ROIs','ROI_coord_ori_lin','ROI_coord_ori_vx_rcp');
    n_ROIs=length(ROI_coord_ori_lin);
    [sz_ori, rez_ori] = read_nrrd_metadata([pathname_raw_anat fish_id '_anatomy.nrrd']);
    affine_x=[pathname_morphed_anat fish_id '_anatomy.xform'];
    ROIs_morphed=zeros(sz_ref,'uint16');
    for i=1:n_ROIs
        % original indeces in voxels
        rcp_ori_vx = double(ROI_coord_ori_vx_rcp{i});
        % original indeces in microns
        rcp_ori_um = rcp_ori_vx.*rez_ori;
        rows_ori_um=rcp_ori_um(:,1);
        columns_ori_um=rcp_ori_um(:,2);
        planes_ori_um=rcp_ori_um(:,3);
        % extend ROIs in z
        plane_ori_um=unique(planes_ori_um);
        ext_planes=plane_ori_um-1:-1:plane_ori_um-rez_ori(3)+1;
        n_ext_planes=length(ext_planes);
        for ii=ext_planes
            planes_ori_um=[planes_ori_um; ones(length(rows_ori_um),1)*ii]; %#ok<AGROW>
        end
        rows_ori_um=repmat(rows_ori_um,n_ext_planes+1,1);
        columns_ori_um=repmat(columns_ori_um,n_ext_planes+1,1);
        % morphed indeces in microns
        rcp_morphed_um=DM_morph_non_ref_coord([rows_ori_um,columns_ori_um,planes_ori_um],affine_x,path_to_lists);
        % morphed indeces in voxels
        rcp_morphed_vx=rcp_morphed_um./rez_ref;
        rcp_morphed_vx=unique([floor(rcp_morphed_vx);ceil(rcp_morphed_vx)],'rows');
        % remove out-of-stack voxels
        bad_coord=find(any(rcp_morphed_vx<=0,2) | rcp_morphed_vx(:,1)>sz_ref(1) | rcp_morphed_vx(:,2)>sz_ref(2) | rcp_morphed_vx(:,3)>sz_ref(3));
        rcp_morphed_vx(bad_coord,:)=[];
        % remove 1/2 voxels that intercept with other already morphed ROIs
        if ~isempty(rcp_morphed_vx)
            % morphed linear coordinates
            lin_morphed=sub2ind(sz_ref, rcp_morphed_vx(:,1), rcp_morphed_vx(:,2), rcp_morphed_vx(:,3));
            % tf stack of this ROI
            this_ROI = false(sz_ref);
            this_ROI(lin_morphed)=true;
            % find if there are ROIs already
            ROIs_there = unique(ROIs_morphed(this_ROI))';
            ROIs_there(ROIs_there==0)=[];
            for ii=ROIs_there
                % find pixels that are intesepting and remove 1/2
                px = find(this_ROI & ROIs_morphed==ii);
                this_ROI(px(randperm(length(px),round(length(px)/2))))=false;
            end
            % add this to ROIs morphed
            ROIs_morphed(this_ROI)=i;
        end
        progressbar([],i/n_ROIs,[]);
    end
    
    % find ROI coordinates
    ROI_coord_morphed_lin=cell(n_ROIs,1);
    for i=1:n_ROIs
        ROI_coord_morphed_lin{i}=find(ROIs_morphed==i);
        progressbar([],[],i/n_ROIs);
    end
    
    % save
    ROI_coord = ROI_coord_morphed_lin;
    save([pathname_ROIs fish_id '_ROIs.mat'],'ROI_coord');
    progressbar(f/n_fish,[],[]);
end