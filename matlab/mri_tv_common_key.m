function key = mri_tv_common_key(cfg)
% Fields that MUST agree when applying frozen calibration parameters.
names={'packageVersion','backend','maxIter','calibrationMaxIter','checkEvery', ...
    'tolKKT','tolGap','tolFeas','noiseStd','clipObservation','calibrationSeeds', ...
    'gridName','diagnosticTau','diagnosticSigma','methodKinds','warmupIterations'};
key=struct();
for j=1:numel(names), key.(names{j})=cfg.(names{j}); end
end
