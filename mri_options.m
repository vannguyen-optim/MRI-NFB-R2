function opts = mri_options(user)
%MRI_OPTIONS Common user settings. Set once and reuse for every stage.
% Examples: opts=mri_options(); opts.outputRoot='D:\Research\MRI_new';
% Changing model, tolerance, or grid requires a new output root and retuning.
root=fileparts(mfilename('fullpath'));
opts.outputRoot=root; % publication: outputs match manuscript paths
opts.gridName='publication48';  % Equal candidate count for all four methods
opts.lambdaTV=3;
opts.noiseStd=10;
opts.maxIter=20000;
opts.tolKKT=1e-4;
opts.tolGap=1e-4;
opts.tolFeas=1e-8;
opts.runSensitivity=true;
opts.compareBackends=true;
opts.makePlots=true;
opts.requireSuccessfulTuning=true;
if nargin<1 || isempty(user), return; end
assert(isstruct(user) && isscalar(user),'MRI:Options','Options must be a scalar struct.');
names=fieldnames(user);
for j=1:numel(names)
    assert(isfield(opts,names{j}),'MRI:Options','Unknown option: %s',names{j});
    opts.(names{j})=user.(names{j});
end
opts.outputRoot=char(opts.outputRoot); opts.gridName=lower(char(opts.gridName));
assert(~isempty(opts.outputRoot),'MRI:Options','Output root cannot be empty.');
assert(any(strcmp(opts.gridName,{'publication48','balanced12','paper12'})),'MRI:Options','Unknown gridName.');
for nm={'lambdaTV','tolKKT','tolGap','tolFeas'}
    validateattributes(opts.(nm{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(opts.maxIter,{'numeric'},{'scalar','real','finite','integer','positive'});
validateattributes(opts.noiseStd,{'numeric'},{'scalar','real','finite','nonnegative'});
for nm={'runSensitivity','compareBackends','makePlots','requireSuccessfulTuning'}
    validateattributes(opts.(nm{1}),{'logical'},{'scalar'});
end
end
