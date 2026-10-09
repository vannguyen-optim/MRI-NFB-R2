function T = run_mri_tuning(imageFiles,options)
%RUN_MRI_TUNING Calibrate the selected fixed grid per method on seeds 1001:1003.
% No test seeds are run. Writes full-precision frozen parameters to tuning.mat.
if nargin<1 || isempty(imageFiles), imageFiles='bird.tif'; end
if nargin<2, options=struct(); end
setup_mri(); opts=mri_options(options); cfg=mri_stage_config('tuning',opts);
T=run_mri_tv(imageFiles,opts.lambdaTV,0,cfg);
fprintf('Frozen parameters: %s\n',fullfile(cfg.outDir,'selected_parameters.csv'));
fprintf('Full calibration audit: %s\n',fullfile(cfg.outDir,'calibration_scores.csv'));
end
