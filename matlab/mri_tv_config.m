function cfg = mri_tv_config(profile)
%MRI_TV_CONFIG Presets for corrected anisotropic-TV experiments.
% Core solver uses only base MATLAB. Pilot/tuning/publication are separate.
if nargin<1, profile='benchmark'; end
cfg.profile=char(profile);
cfg.packageVersion='MRI-TV-2026-10-08-publication48';
cfg.outDir='';
cfg.backend='spatial';
cfg.maxIter=20000;
cfg.calibrationMaxIter=20000;
cfg.checkEvery=10;
cfg.tolKKT=1e-4;
cfg.tolGap=1e-4;
cfg.tolFeas=1e-8;
cfg.noiseStd=10;
cfg.clipObservation=true;
cfg.seedBase=2026;
cfg.calibrationSeeds=[1001 1002 1003];
cfg.gridName='publication48';
cfg.tune=true;
cfg.tuningOnly=false;
cfg.tuningFile='';
cfg.requireSuccessfulTuning=true;
cfg.sensitivity=false;
cfg.sensitivityAlpha=[0 .01 .05 .1 .2 .25 .3 .6 .9 1.2 2];
cfg.sensitivityRho=[.1 .2 .3];
cfg.sensitivityTau=.2;
cfg.sensitivitySigma=.3;
cfg.compareBackends=false;
cfg.makePlots=true;
cfg.saveImages=true;
cfg.warmupIterations=5;
cfg.runTests=true;
cfg.resume=true;
cfg.diagnosticTau=.2;
cfg.diagnosticSigma=.3;
cfg.methodKinds={'MRI','RI','CV','CP'};
switch lower(char(profile))
    case 'smoke'
        cfg.maxIter=100; cfg.calibrationMaxIter=100;
        cfg.tune=false; cfg.seedBase=8000; cfg.makePlots=false;
    case 'pilot'
        cfg.tune=false; cfg.seedBase=9000;
    case 'tuning'
        cfg.tuningOnly=true;
    case 'publication'
        cfg.sensitivity=true; cfg.compareBackends=true;
    case 'benchmark'
        % Low-level backward-compatible profile, with 3 calibration seeds.
    otherwise
        error('MRI:Profile','Use smoke, pilot, tuning, publication, or benchmark.');
end
end
