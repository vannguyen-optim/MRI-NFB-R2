function imageFile = mri_select_image(imageFile)
%MRI_SELECT_IMAGE Resolve the reference image without changing its pixels.
root=setup_mri();
if nargin<1,imageFile='';end
if isempty(imageFile)
    candidates={fullfile(root,'bird.tif'),fullfile(root,'data','bird.tif')};
    for j=1:numel(candidates)
        if exist(candidates{j},'file')==2,imageFile=candidates{j};break;end
    end
end
if isempty(imageFile)
    if ~usejava('desktop')
        error('MRI:MissingImage','Put bird.tif in the package root/data folder or supply its full path.');
    end
    [name,folder]=uigetfile( ...
        {'*.tif;*.tiff;*.png;*.bmp','Grayscale reference image'}, ...
        'Select the reference image for Experiment 2',root);
    if isequal(name,0),error('MRI:MissingImage','Image selection was cancelled.');end
    imageFile=fullfile(folder,name);
end
imageFile=char(imageFile);
[ok,info]=fileattrib(imageFile);
assert(ok && ~info.directory,'MRI:MissingImage','Image not found: %s',imageFile);
imageFile=info.Name;
end
