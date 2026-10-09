%RUN_EXPERIMENT_2 Subsection 4.2: TV IMAGE RESTORATION only.
% Includes pilot, tuning, 20 tests, sensitivity, direct/FFT and report export.
% Put bird.tif in this folder or data/, or choose it in the dialog.
root=fileparts(mfilename('fullpath'));
cd(root);addpath(root,'-begin');
EXPERIMENT2_RESULTS=run_tv_publication_all(); 
