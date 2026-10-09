function report = export_mri_results(resultsDir,makePlots)
%EXPORT_MRI_RESULTS Re-export tables and figures without rerunning solvers.
% Point at a directory containing benchmark.mat, not at tuning/.
root=setup_mri();
if nargin<1 || isempty(resultsDir), resultsDir=fullfile(root,'results_tv'); end
if nargin<2, makePlots=true; end
report=export_mri_tv_paper(char(resultsDir),makePlots);
end
