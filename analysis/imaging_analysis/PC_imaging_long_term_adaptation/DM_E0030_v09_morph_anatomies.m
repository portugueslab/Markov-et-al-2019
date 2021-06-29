% this code is not supposed to be executed as a whole. Run sections
% separately and check manually how thigs go

%% initial stuff
clc; close all; clear all;
pathname_results='C:\Users\dmarkov\Desktop\Data_Daniil\E0030_v09\';
pathname_raw_anat=[pathname_results 'anatomies\raw_nrrds\'];
pathname_final_anat = [pathname_results 'anatomies\morphed\'];
pathname_MS_data = '...\Markov_et_al_2021_Nat_Commun_data&code\data\';
pathname_ref = [pathname_MS_data 'reference_brain_stacks\'];
filename_ref = 'PortuguesLab_PC_ref.nrrd';
str_ref=[pathname_ref filename_ref];
sz_ref = read_nrrd_metadata(str_ref);

%% try to morph all files to the original reference
cd(pathname_raw_anat);
all_fish=dm_dir('*_f*_anatomy.nrrd');
progressbar('Morphing anatomies...');
for f=1:length(all_fish)
    str_final = [pathname_final_anat all_fish{f}];
    if exist(str_final,'file')~=2
        str_morph = [pathname_raw_anat all_fish{f}];
        affine_x = strrep(str_final,'.nrrd','.xform');
        morph_stack_to_ref(str_ref, str_morph, str_final, affine_x);
    end
    progressbar(f/length(all_fish));
end
% this resulted in 15 acceptable morphings. Manually delete 10 bad ones and
% compute the mean transformed stack

%% create average stack
cd(pathname_final_anat);
all_fish=dm_dir('*_f*_anatomy.nrrd');
A=zeros(sz_ref);
for f=1:length(all_fish)
    str_final = [pathname_final_anat all_fish{f}];
    B = double(nrrdread(str_final));
    A=A+B;
end
A=uint8(A/max(A(:))*255);
mat2tiff(A,pathname_final_anat,['average_' num2str(length(all_fish)) '_brains.tif']);
% set the correct rezolution manually and save as nrrd

%% now try to morph the rest of the fish to this new averaged stack
str_ref2=[pathname_final_anat 'average_24_brains.nrrd'];
cd(pathname_raw_anat);
all_fish=dm_dir('*_f*_anatomy.nrrd');
progressbar('Morphing anatomies...');
for f=1:length(all_fish)
    str_final = [pathname_final_anat all_fish{f}];
    if exist(str_final,'file')~=2
        str_morph = [pathname_raw_anat all_fish{f}];
        affine_x = strrep(str_final,'.nrrd','.xform');
        morph_stack_to_ref(str_ref2, str_morph, str_final, affine_x);
    end
    progressbar(f/length(all_fish));
end
% morphing to average of 15 brains: one more fish. Create another averaged stack and try again
% morphing to average of 19 brains: no new fish. Change the approach
% morphing to average of 24 brains: did not work for the last stack

%% Try to morph each of the remaining fish to each of the morphed fish
cd(pathname_raw_anat);
all_fish=dm_dir('*_f*_anatomy.nrrd');
all_morphed_fish=dm_dir([pathname_final_anat '*_f*_anatomy.nrrd']);
progressbar('Morphing anatomies...');
for f=1:length(all_fish)
    if exist([pathname_final_anat all_fish{f}],'file')~=2
        str_morph = [pathname_raw_anat all_fish{f}];
        pathname_final2=[pathname_final_anat strrep(all_fish{f},'_anatomy.nrrd','') '_individual_morphings\'];
        mkdir(pathname_final2);
        for j=1:length(all_morphed_fish)
            str_ref2 = [pathname_final_anat all_morphed_fish{j}];
            str_final = [pathname_final2 strrep(all_fish{f},'_anatomy.nrrd','') '_2_' strrep(all_morphed_fish{j},'_anatomy.nrrd','') '.nrrd'];
            affine_x = strrep(str_final,'.nrrd','.xform');
            morph_stack_to_ref(str_ref2, str_morph, str_final, affine_x);
        end
    end
    progressbar(f/length(all_fish));
end
% morph 9 fish to 16 morphed fish: resulted in OK morphings for 8 fish. Keep the best nrrd stack
% morph the remaining fish to 24 morphed fish: kinda worked

%% move them to the original folder
cd(pathname_final_anat);
all_fish=dm_dir('*_f*_individual_morphings');
for f=1:length(all_fish)
    fish_id=strrep(all_fish{f},'_individual_morphings','');
    cd([pathname_final_anat all_fish{f}]);
    ref_fish_id = dm_dir('*_f*_2_*_f*.nrrd');
    movefile([pathname_final_anat fish_id '_individual_morphings\' ref_fish_id{1}], [pathname_final_anat fish_id '_anatomy.nrrd']);
    movefile(strrep([pathname_final_anat fish_id '_individual_morphings\' ref_fish_id{1}],'.nrrd','.xform'), [pathname_final_anat fish_id '_anatomy.xform']);
end

%% to check the everything is correct, apply computed transformations to each stack
cd(pathname_raw_anat);
all_fish=dm_dir('*_f*_anatomy.nrrd');
progressbar('Morphing anatomies...');
for f=1:length(all_fish)
    str_final = [pathname_final_anat all_fish{f}];
    if exist(str_final,'file')~=2
        str_morph = [pathname_raw_anat all_fish{f}];
        affine_x = strrep(str_final,'.nrrd','.xform');
        apply_morphing(str_ref, str_morph, str_final,affine_x)
    end
    progressbar(f/length(all_fish));
end

% final 25 morphed anatomies are the ones uploaded with this manuscript into this folder:
% ...\Markov_et_al_2021_Nat_Commun_data&code\data\PC_imaging_long_term_adaptation\anatomies