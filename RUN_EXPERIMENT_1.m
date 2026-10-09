%RUN_EXPERIMENT_1 Subsection 4.1: SYNTHETIC nonlinear varying-kernel example.
root=fileparts(mfilename('fullpath'));
cd(root);addpath(root,'-begin');
EXPERIMENT1_RESULTS=run_mri_synthetic();
