function report = mri_preflight(imageFiles,options,verbose)
%MRI_PREFLIGHT Check image, file access and all predetermined candidates.
% No solver is run and no experiment result is created by this check.
setup_mri();
if nargin<3, verbose=true; end
if nargin<2, options=struct(); end
options=mri_options(options);
[truth,names]=mri_tv_load_images(imageFiles);
folder=options.outputRoot;
if exist(folder,'dir')~=7
    [ok,msg]=mkdir(folder);
    assert(ok,'MRI:Output','Cannot create output directory: %s',msg);
end
probe=[tempname(folder) '.tmp'];
fid=fopen(probe,'w');
assert(fid>=0,'MRI:Output','Cannot write to the output directory: %s',folder);
fprintf(fid,'MRI-NFB write-access test\n'); fclose(fid); delete(probe);
% Mathematical candidate checks are independent of the image or test data.
methods={'MRI','RI','CV','CP'}; counts=zeros(1,numel(methods));
for j=1:numel(methods)
    candidates=mri_tv_parameters(methods{j},true,options.gridName);
    counts(j)=numel(candidates);
    assert(counts(j)>0,'MRI:Grid','The candidate grid is empty.');
    for k=1:numel(candidates)
        mri_tv_check_parameters(candidates(k));
    end
end
assert(all(counts==counts(1)),'MRI:Grid','Methods must have equal candidate counts.');
report=struct('passed',true,'imageNames',{names},'imageSizes', ...
    {cellfun(@size,truth,'UniformOutput',false)},'outputRoot',folder, ...
    'gridName',options.gridName,'methods',{methods},'candidateCounts',counts);
if ~verbose, return; end
fprintf('\nMRI-NFB run configuration\n');
for j=1:numel(truth)
    fprintf('  Image %d: %s (%d x %d)\n',j,names{j},size(truth{j},1),size(truth{j},2));
end
fprintf('  Output: %s\n',folder);
fprintf('  Grid: %s | %d candidates/method | fixed lambdaTV = %.8g\n', ...
    options.gridName,counts(1),options.lambdaTV);
fprintf('  Pilot seeds: 9001:9003 | tuning seeds: 1001:1003 | test seeds: 2027:2046\n');
fprintf('  MaxIter = %d | KKT = %.3g | gap = %.3g | feasibility = %.3g\n', ...
    options.maxIter,options.tolKKT,options.tolGap,options.tolFeas);
end
