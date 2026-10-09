function OUT = run_publication_all(imageFile)
%RUN_PUBLICATION_ALL Both experiments in the manuscript, executed once each.
% 1: nonlinear varying-kernel synthetic example (no image needed).
% 2: TV restoration with independent pilot, tuning, tests and report export.
if nargin<1,imageFile='';end
root=setup_mri();opts=mri_options();
imageFile=mri_select_image(imageFile);
mri_preflight(imageFile,opts,true);
checks=check_mri_installation(true);
assert(checks.passed,'MRI:Checks','The installation checks did not pass.');
fprintf('\nEXPERIMENT 1 / 2: synthetic nonlinear varying kernels\n');
S=run_mri_synthetic(opts);
assert(all(S.summary.Converged),'MRI:SyntheticCap','Inspect the synthetic summary.');
assert(all(S.summary.MaxBackwardResidual<=1e-13), ...
    'MRI:SyntheticBackward','A backward solve exceeded its prescribed tolerance.');
fprintf('\nEXPERIMENT 2 / 2: TV-regularized image restoration\n');
TV=run_tv_publication_all(imageFile,opts,false);
OUT=struct('checks',checks,'synthetic',S,'tv',TV,'options',opts);
save(fullfile(root,'publication_run_summary.mat'),'OUT','-v7');
fprintf('\nBoth experiments completed.\n');
fprintf('Experiment 1: %s\n',fullfile(opts.outputRoot,'results_synthetic'));
fprintf('Experiment 2: %s\n',fullfile(opts.outputRoot,'results_tv'));
end
