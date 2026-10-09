function [truth,names] = mri_tv_load_images(imageFiles)
% Read without silent resizing, RGB conversion, rescaling, or replacement.
if isnumeric(imageFiles), sources={imageFiles};
elseif ischar(imageFiles), sources={imageFiles};
elseif isstring(imageFiles), sources=cellstr(imageFiles);
elseif iscell(imageFiles), sources=imageFiles;
else, error('MRI:Images','Use a grayscale array, filename, or cell of filenames.'); end
assert(~isempty(sources),'MRI:Images','At least one image is required.');
truth=cell(1,numel(sources)); names=cell(size(truth));
for i=1:numel(sources)
    if isnumeric(sources{i})
        x=double(sources{i}); names{i}=sprintf('array_%d',i);
    else
        names{i}=char(sources{i});
        if ~exist(names{i},'file'), error('MRI:MissingImage','Image not found: %s',names{i}); end
        [x,map]=imread(names{i});
        assert(isempty(map),'MRI:ImageFormat','Use grayscale, not indexed, images.');
        assert(ndims(x)==2,'MRI:ImageFormat','Supply a grayscale image; RGB is not converted silently.');
        assert(isa(x,'uint8') || isfloat(x),'MRI:ImageScale', ...
            'Use uint8 or explicit floating-point intensities in [0,255].');
        x=double(x);
    end
    validateattributes(x,{'double'},{'2d','real','finite','nonempty','>=',0,'<=',255});
    assert(min(size(x))>=2,'MRI:ImageSize','Each dimension must be at least 2.');
    assert(max(x(:))>1,'MRI:ImageScale','Use [0,255] intensities, not [0,1].');
    truth{i}=x;
end
end
