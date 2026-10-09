%RUN_PUBLICATION Run BOTH manuscript experiments, using MATLAB only.
% Experiment 1: synthetic. Experiment 2: TV restoration.
% Place bird.tif in the package root or data/; no image is resized or replaced.
root=fileparts(mfilename('fullpath'));
cd(root);addpath(root,'-begin');
PUBLICATION_RESULTS=run_publication_all(); %#ok<NASGU>
