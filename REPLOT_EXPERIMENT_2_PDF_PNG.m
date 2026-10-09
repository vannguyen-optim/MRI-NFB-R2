function report = REPLOT_EXPERIMENT_2_PDF_PNG(resultsDir)
%REPLOT_EXPERIMENT_2_PDF_PNG Export saved TV figures without running solvers.
% Place this file next to RUN_EXPERIMENT_2.m. Replace the supplied file
% reporting/export_mri_tv_paper.m before using this replot command.
%
% REPLOT_EXPERIMENT_2_PDF_PNG
% REPLOT_EXPERIMENT_2_PDF_PNG('full/path/to/results_tv')
%
% Reads benchmark.mat. Writes only figure files in results_tv/figures/.
% Does not alter CSV/MAT data, reselect parameters, or rerun any solver.
root=fileparts(mfilename('fullpath'));
assert(exist(fullfile(root,'setup_mri.m'),'file')==2,'MRI:PlotInstall', ...
    'Place this file in the package root next to RUN_EXPERIMENT_2.m.');
addpath(root,'-begin');
setup_mri();
if nargin<1 || isempty(resultsDir)
    resultsDir=fullfile(root,'results_tv');
    if exist(fullfile(resultsDir,'benchmark.mat'),'file')~=2
        alternate=fullfile(root,'results','results_tv');
        if exist(fullfile(alternate,'benchmark.mat'),'file')==2
            resultsDir=alternate;
        else
            [name,folder]=uigetfile({'benchmark.mat','Completed TV results (benchmark.mat)'}, ...
                'Select benchmark.mat from results_tv',root);
            if isequal(name,0)
                report=struct('cancelled',true,'solverRerun',false);
                return;
            end
            assert(strcmpi(name,'benchmark.mat'),'MRI:PlotResults', ...
                'Select the benchmark.mat file from the completed TV experiment.');
            resultsDir=folder;
        end
    end
end
resultsDir=char(resultsDir);
if exist(resultsDir,'file')==2
    [folder,name,extension]=fileparts(resultsDir);
    assert(strcmpi([name extension],'benchmark.mat'),'MRI:PlotResults', ...
        'Use results_tv/ or its benchmark.mat file.');
    resultsDir=folder;
    if isempty(resultsDir), resultsDir=pwd; end
end
assert(exist(fullfile(resultsDir,'benchmark.mat'),'file')==2,'MRI:PlotResults', ...
    'No benchmark.mat found in %s.',resultsDir);
assert(nargin('export_mri_tv_paper')==3,'MRI:PlotVersion', ...
    'Replace reporting/export_mri_tv_paper.m with the accompanying version.');
report=export_mri_tv_paper(resultsDir,true,true);
end
