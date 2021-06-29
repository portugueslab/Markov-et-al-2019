function [] = mat2tiff(mat_file, pathname, filename)
warning('off','MATLAB:DELETE:FileNotFound');
delete(fullfile(pathname, filename))
if length(size(mat_file))==3
    num_planes=size(mat_file,3);
    progressbar('saving tif stack')
    for plane=1:num_planes
        imwrite(mat_file(:,:,plane),fullfile(pathname, filename),'writemode','append');
        progressbar(plane/num_planes);
    end
elseif length(size(mat_file))==4
    sz=size(mat_file);
    if sz(3)~=3
        error('3d dimension of the stack should be rgb, not planes!');
    else
        num_planes=size(mat_file,4);
        progressbar('saving tif stack')
        for plane=1:num_planes
            imwrite(mat_file(:,:,:,plane),fullfile(pathname, filename),'writemode','append');
            progressbar(plane/num_planes);
        end
    end
else
    error('Unexpected number of dimensions in the stack. It should be 3 or 4');
end