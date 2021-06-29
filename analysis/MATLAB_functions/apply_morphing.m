function [] = apply_morphing(str_ref, str_morph, str_final,affine_x)

command_str_reform=['C:\cmtk_files\bin\reformatx --echo -o' ' ' str_final ' ' '--floating' ' ' str_morph ' ' str_ref ' ' affine_x];
dos(command_str_reform,'-echo');