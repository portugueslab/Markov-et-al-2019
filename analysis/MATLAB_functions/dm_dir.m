function [cell_output] = dm_dir(str)
cell_output=struct2cell(dir(str));
cell_output=cell_output(1,:)';