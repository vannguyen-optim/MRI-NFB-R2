function P = run_mri_pilot(imageFiles,options)
%RUN_MRI_PILOT Three independent pilot seeds 9001:9003, without tuning.
% Uses the full supplied image (no resizing) and the publication tolerances.
% P=run_mri_pilot('bird.tif');
if nargin<1 || isempty(imageFiles), imageFiles='bird.tif'; end
if nargin<2, options=struct(); end
setup_mri(); opts=mri_options(options); cfg=mri_stage_config('pilot',opts);
P=run_mri_tv(imageFiles,opts.lambdaTV,3,cfg);
export_mri_results(cfg.outDir,opts.makePlots);
for kind=cfg.methodKinds
    rr=P.records(cellfun(@(z)strcmp(z.method.kind,kind{1}),P.records));
    fprintf('PILOT %s: %d/%d successful.\n',kind{1},sum(cellfun(@(z)z.success,rr)),numel(rr));
end
end
