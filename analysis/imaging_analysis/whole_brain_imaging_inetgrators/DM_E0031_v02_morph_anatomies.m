% this will morph anatomies to the reference brain
clc; clear all; close all;
addpath('J:\Daniil Markov\matlab\real_functions\morphing');
pathname_morph='C:\Users\dmarkov\Desktop\Data_Daniil\Internal_models_figures\morphing\flipped_anatomies\';
pathname_final='C:\Users\dmarkov\Desktop\Data_Daniil\Internal_models_figures\morphing\flipped_anatomies_morphed\';
pathname_ref='C:\Users\dmarkov\Desktop\Data_Daniil\Internal_models_figures\morphing\';
filename_ref='lightsheet_ref_morphed_to_PortuguesLab_wholebrain_ref';
str_ref=[pathname_ref filename_ref];

%% morph all stacks to ref
cd(pathname_morph);
all_fish=dm_dir('*_f*_anatomy.nrrd');
progressbar('Morphing anatomies...')
for f=1:length(all_fish)
    this_name=all_fish{f};
    str_final=[pathname_final this_name];
    if exist(str_final,'file')~=2
        str_morph=[pathname_morph this_name];
        affine_x=strrep(str_final,'.nrrd','.xform');
        morph_stack_to_ref(str_ref, str_morph, str_final,affine_x);
    end
    progressbar(f/length(all_fish));
end
% this resulted in 3 good fish, manualy delete the other 3 from the folder

%% now let's try to morph the other 2 fish
all_fish ={'190524_f0_anatomy.nrrd';'190524_f1_anatomy.nrrd';'190524_f2_anatomy.nrrd'};

%% things that didn't work:
% 1. morphing to the same ref with cropping boxes
% 2. morphing to average 4 good fish with cropping boxes
% 3. morphing to each morphed good fish
% 4. morphing with aninitial transformation taken from random good fish
% 5. mrphing to each raw good fish
% 6. changing initial rezolution from 0.6 to 0.65
% 7. morphing to ls fish that was also imaged under the confocal
% 8. morphing using initial transformation

%% try to  manually edit the xfrom folder and see what happens
% first, just morph as usual to create some xfrom folders
for f=1:length(all_fish)
    this_name=all_fish{f};
    str_final=[pathname_final this_name];
    str_morph=[pathname_morph this_name];
    affine_x_init=strrep(str_final,'.nrrd','.xform');
    morph_stack_to_ref(str_ref, str_morph, str_final,affine_x_init);
end

%% now manually play with registration parameters for both fish
this_name=all_fish{3};
str_final=[pathname_final strrep(this_name,'anatomy','anatomy_init_tranform')];
str_morph=[pathname_morph this_name];
affine_x_init=[pathname_final strrep(this_name,'anatomy.nrrd','anatomy_init.xform')];
apply_morphing(str_ref, str_morph, str_final, affine_x_init);

%% morph them to the reference using manually created initial transformation
this_name=all_fish{3};
str_final=[pathname_final this_name];
str_morph=[pathname_morph this_name];
affine_x=strrep(str_final,'.nrrd','.xform');
affine_x_init=[pathname_final strrep(this_name,'anatomy.nrrd','anatomy_init.xform')];
morph_stack_to_ref(str_ref, str_morph, str_final,affine_x,[],[],affine_x_init);

%% now make sure that the xfrom folders are correct
cd(pathname_morph);
all_fish=dm_dir('*_f*_anatomy.nrrd');
progressbar('Morphing anatomies...')
for f=1:length(all_fish)
    this_name=all_fish{f};
    str_final=[pathname_final this_name];
    str_morph=[pathname_morph this_name];
    affine_x=strrep(str_final,'.nrrd','.xform');
    apply_morphing(str_ref, str_morph, str_final,affine_x);
    progressbar(f/length(all_fish));
end

% final 6 morphed anatomies are the ones uploaded with this manuscript into this folder:
% ...\Markov_et_al_2021_Nat_Commun_data&code\data\whole_brain_imaging_inetgrators\anatomies