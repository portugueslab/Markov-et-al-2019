function [] = morph_stack_to_ref(str_ref, str_morph, str_final, affine_x, ref_cropping_box, float_cropping_box, affine_x_init)
% cropping_box: 'x0,y0,z0,x1,y1,z1' in pixels
path_to_cmtk='C:\cmtk_files\bin\';
if nargin == 4
    command_str_affine=[path_to_cmtk 'registration --echo --initxlate --dofs 6,9 --auto-multi-levels 4 --match-histograms --exploration 150 --accuracy 0.1 --coarsest 8 --omit-original-data -o' ' ' affine_x ' ' str_ref ' ' str_morph]; 
elseif nargin == 5
    command_str_affine=[path_to_cmtk 'registration --echo --crop-index-ref ' ref_cropping_box ' --initxlate --dofs 6,9 --auto-multi-levels 4 --match-histograms --exploration 150 --accuracy 0.1 --coarsest 8 --omit-original-data -o' ' ' affine_x ' ' str_ref ' ' str_morph];
elseif nargin == 6
    command_str_affine=[path_to_cmtk 'registration --echo --crop-index-ref ' ref_cropping_box ' --crop-index-flt ' float_cropping_box ' --initxlate --dofs 6,9 --auto-multi-levels 4 --match-histograms --exploration 150 --accuracy 0.1 --coarsest 8 --omit-original-data -o' ' ' affine_x ' ' str_ref ' ' str_morph];
else
    command_str_affine=[path_to_cmtk 'registration --echo --crop-index-ref ' ref_cropping_box ' --crop-index-flt ' float_cropping_box ' --initial ' affine_x_init ' --initxlate --dofs 6,9 --auto-multi-levels 4 --match-histograms --exploration 150 --accuracy 0.1 --coarsest 8 --omit-original-data -o' ' ' affine_x ' ' str_ref ' ' str_morph];
    if isempty(ref_cropping_box) && isempty(float_cropping_box)
        command_str_affine=[path_to_cmtk 'registration --echo --initial ' affine_x_init ' --initxlate --dofs 6,9 --auto-multi-levels 4 --match-histograms --exploration 150 --accuracy 0.1 --coarsest 8 --omit-original-data -o' ' ' affine_x ' ' str_ref ' ' str_morph];
    else
        error('fix it!');
    end
end
command_str_reform=[path_to_cmtk 'reformatx --echo -o' ' ' str_final ' ' '--floating' ' ' str_morph ' ' str_ref ' ' affine_x];
dos(command_str_affine,'-echo');
dos(command_str_reform,'-echo');