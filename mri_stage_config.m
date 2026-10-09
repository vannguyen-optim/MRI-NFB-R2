function cfg = mri_stage_config(stage,opts)
%MRI_STAGE_CONFIG No decisions are made using test-set results.
cfg=mri_tv_config(stage);
cfg.outDir=fullfile(opts.outputRoot,stage);
cfg.gridName=opts.gridName; cfg.noiseStd=opts.noiseStd;
if ~strcmp(stage,'smoke')
    cfg.maxIter=opts.maxIter; cfg.calibrationMaxIter=opts.maxIter;
end
cfg.tolKKT=opts.tolKKT; cfg.tolGap=opts.tolGap; cfg.tolFeas=opts.tolFeas;
cfg.makePlots=opts.makePlots;
cfg.requireSuccessfulTuning=opts.requireSuccessfulTuning;
if strcmp(stage,'publication')
    cfg.outDir=fullfile(opts.outputRoot,'results_tv');
    cfg.tuningFile=fullfile(opts.outputRoot,'tuning','tuning.mat');
    cfg.sensitivity=opts.runSensitivity;
    cfg.compareBackends=opts.compareBackends;
end
end
