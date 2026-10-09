function R = run_mri_publication(imageFiles,options)
%RUN_MRI_PUBLICATION Run the frozen-parameter TV publication stage.
% This function does NOT retune after seeing test results and does not rerun
% the synthetic experiment. If tuning.mat is absent, calibration is performed
% first on the disjoint seeds 1001:1003.
%
% R = run_mri_publication('bird.tif');
%
% Publication test seeds are 2027:2046. The stage also runs the sensitivity
% panel and direct/FFT implementation comparison, then exports the paper
% tables and figures.
if nargin<1 || isempty(imageFiles), imageFiles='bird.tif'; end
if nargin<2, options=struct(); end
setup_mri();
opts=mri_options(options);
cfg=mri_stage_config('publication',opts);
if ~exist(cfg.tuningFile,'file')
    T=run_mri_tuning(imageFiles,opts); %#ok<NASGU>
end
R=run_mri_tv(imageFiles,opts.lambdaTV,20,cfg);
export_mri_results(cfg.outDir,opts.makePlots);
end
