function [ mat_file ] = tiff2mat( pathname, filename )
%this function transorms a tiff stack in a file to a mat file in the workspace
num_planes=length(imfinfo(fullfile(pathname, filename)));
test_plane=imread(fullfile(pathname, filename),1);
mat_file=zeros([size(test_plane) num_planes],  class(test_plane));
if ndims(mat_file)==3 && num_planes>1
    for plane=1:num_planes
        mat_file(:,:,plane)=imread(fullfile(pathname, filename),plane);
    end
elseif ndims(mat_file)==4 || (ndims(mat_file)==3 && num_planes==1)
    for plane=1:num_planes
        mat_file(:,:,:,plane)=imread(fullfile(pathname, filename),plane);
    end
else
    error('Strange number of dimensions');
end

