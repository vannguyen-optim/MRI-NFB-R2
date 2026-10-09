function OUT = run_tv_publication_all(imageFile,options,runChecks)
%RUN_TV_PUBLICATION_ALL EXPERIMENT 2 only: pilot, tuning, test and export.
% No synthetic run is performed here. No Python is required.
if nargin<1,imageFile='';end
if nargin<2,options=struct();end
if nargin<3,runChecks=true;end
setup_mri();opts=mri_options(options);imageFile=mri_select_image(imageFile);
mri_preflight(imageFile,opts,true);
if runChecks
    checks=check_mri_installation(true);
    assert(checks.passed,'MRI:Checks','The installation checks did not pass.');
end
fprintf('\nEXPERIMENT 2: TV-restoration pilot\n');
P=run_mri_pilot(imageFile,opts);
assert(~any(cellfun(@(r)strcmp(r.status,'numerical_failure'),P.records)), ...
    'MRI:PilotFailure','Inspect pilot/raw_results.csv: a numerical failure occurred.');
fprintf('\nEXPERIMENT 2: calibration on seeds 1001:1003\n');
T=run_mri_tuning(imageFile,opts);
assert(T.allCalibrationSuccessful,'MRI:CalibrationFailed', ...
    'Selected parameters must pass all three calibration runs. See calibration_scores.csv.');
fprintf('\nEXPERIMENT 2: 20 frozen-parameter tests, sensitivity and direct/FFT\n');
R=run_mri_publication(imageFile,opts);
OUT=struct('pilot',P,'tuning',T,'tv',R,'imageFile',imageFile,'options',opts);
save(fullfile(opts.outputRoot,'experiment2_summary.mat'),'OUT','-v7');
fprintf('\nExperiment 2 completed: %s\n',fullfile(opts.outputRoot,'results_tv'));
end
