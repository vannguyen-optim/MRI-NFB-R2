function R = run_mri_tv(imageFiles,lambdaTV,nRuns,cfg)
%RUN_MRI_TV Corrected TV experiment with separate calibration and test stages.
% Use root-level run_mri_pilot, run_mri_tuning, run_mri_publication for normal use.
% All four methods solve the SAME objective; no legacy numerical data enter.
if nargin<2 || isempty(lambdaTV), error('MRI:Weight','Supply lambdaTV (reference value: 3).'); end
if nargin<3 || isempty(nRuns), nRuns=20; end
if nargin<4 || isempty(cfg), cfg=mri_tv_config('benchmark'); end
validate_cfg(cfg,nRuns,lambdaTV);
[truth,names]=mri_tv_load_images(imageFiles); nImages=numel(truth);
if isempty(cfg.outDir), cfg.outDir=fullfile(pwd,['mri_tv_' cfg.profile]); end
cfg.outDir=char(cfg.outDir);
if ~exist(cfg.outDir,'dir'), mkdir(cfg.outDir); end
lock=fullfile(cfg.outDir,'.running_lock');
if exist(lock,'dir'), error('MRI:Locked','Another run or stale lock exists: %s',lock); end
[ok,msg]=mkdir(lock);
if ~ok || ~isempty(msg), error('MRI:Locked','Cannot acquire run lock: %s',lock); end
lockGuard=onCleanup(@()release_lock(lock)); %#ok<NASGU>
for nm={'tasks','tables','figures','data','data_calibration'}
    p=fullfile(cfg.outDir,nm{1}); if ~exist(p,'dir'), mkdir(p); end
end
oldrng=rng; rngGuard=onCleanup(@()rng(oldrng)); %#ok<NASGU>
code=mri_tv_source_snapshot();
[fixed,importedCal,importedScores,tuningInfo]=read_frozen(cfg,truth,lambdaTV,code);
key=rmfield(cfg,{'outDir','resume','runTests','makePlots','saveImages','tuningFile'});
manifestFile=fullfile(cfg.outDir,'manifest.mat');
if exist(manifestFile,'file')
    z=load(manifestFile,'manifest'); manifest=z.manifest;
    assert(cfg.resume,'MRI:Resume','Output exists; resume is disabled. Use a new output directory.');
    assert(isequaln(manifest.cfgKey,key) && isequaln(manifest.truth,truth) && ...
        manifest.lambdaTV==lambdaTV && manifest.nRuns==nRuns, ...
        'MRI:Resume','Configuration/image/weight/run count changed. Use a new output directory.');
    assert(isequaln(manifest.codeSnapshot,code),'MRI:Resume', ...
        'Core source changed. Do not mix results; use a new output directory.');
    assert(strcmp(manifest.version,version) && strcmp(manifest.computer,computer), ...
        'MRI:Resume','MATLAB/platform changed. Use a new output directory.');
    assert(isequaln(manifest.frozenSelection,fixed),'MRI:Resume', ...
        'Frozen parameters changed. Use a new output directory.');
else
    manifest.cfgKey=key; manifest.commonKey=mri_tv_common_key(cfg);
    manifest.truth=truth; manifest.imageNames=names; manifest.lambdaTV=lambdaTV;
    manifest.nRuns=nRuns; manifest.version=version; manifest.computer=computer;
    manifest.codeSnapshot=code; manifest.created=datestr(now,30);
    manifest.frozenSelection=fixed; manifest.tuningInfo=tuningInfo;
    manifest.timing='Sequential tic/toc wall-clock time INCLUDING stopping diagnostics';
    manifest.note='Corrected exact TV model; no legacy table values imported.';
    mri_tv_atomic_save(manifestFile,'manifest',manifest);
end
if cfg.runTests
    checks=test_mri_tv(false);
    mri_tv_atomic_save(fullfile(cfg.outDir,'selftests.mat'),'checks',checks);
end
write_grid(cfg);
records={}; calibration={}; scores=cell(nImages,numel(cfg.methodKinds));
selected=cell(nImages,numel(cfg.methodKinds));
if ~isempty(fixed)
    selected=fixed; calibration=importedCal; scores=importedScores;
end
for im=1:nImages
    xtrue=truth{im}; op=mri_tv_operators(size(xtrue),cfg.backend);
    opData=mri_tv_operators(size(xtrue),'spatial'); clean=opData.K(xtrue);
    if isempty(fixed)
        if cfg.tune
            [chosen,cal,sc]=calibrate_image(im,xtrue,clean,lambdaTV,op,cfg);
            selected(im,:)=chosen; scores(im,:)=sc;
            calibration=[calibration cal]; %#ok<AGROW>
        else
            for j=1:numel(cfg.methodKinds)
                selected{im,j}=mri_tv_parameters(cfg.methodKinds{j},false,cfg.gridName);
            end
        end
    end
    % Write selected/calibration output even for tuning-only or interrupted jobs.
    R=pack(records,selected,calibration,scores,manifest,cfg);
    mri_tv_atomic_save(fullfile(cfg.outDir,'progress.mat'),'R',R);
    mri_tv_export(records,selected,calibration,manifest,cfg,false);
    if cfg.tuningOnly, continue; end
    for r=1:nRuns
        seed=cfg.seedBase+r;
        b=save_observation(clean,seed,cfg,im,false);
        if r==1 && cfg.saveImages
            imwrite(uint8(xtrue),fullfile(cfg.outDir,'figures',sprintf('i%02d_truth.png',im)));
            imwrite(uint8(b),fullfile(cfg.outDir,'figures',sprintf('i%02d_observed.png',im)));
        end
        rng(55000+1000*im+r,'twister'); order=randperm(numel(cfg.methodKinds));
        for ii=order
            m=selected{im,ii}; taskid=sprintf('test_i%02d_%s_s%d',im,m.kind,seed);
            rec=task(taskid,im,seed,b,xtrue,lambdaTV,op,m,cfg,'benchmark',r==1);
            records{end+1}=rec; %#ok<AGROW>
        end
        R=pack(records,selected,calibration,scores,manifest,cfg);
        mri_tv_atomic_save(fullfile(cfg.outDir,'progress.mat'),'R',R);
        mri_tv_export(records,selected,calibration,manifest,cfg,false);
    end
    if cfg.sensitivity
        grid=sensitivity_grid(cfg);
        for r=1:nRuns
            seed=cfg.seedBase+r; b=save_observation(clean,seed,cfg,im,false);
            rng(65000+1000*im+r,'twister'); order=randperm(numel(grid));
            for ii=order
                m=grid{ii}; taskid=sprintf('sens_i%02d_m%03d_s%d',im,m.id,seed);
                rec=task(taskid,im,seed,b,xtrue,lambdaTV,op,m,cfg,'sensitivity',false);
                records{end+1}=rec; %#ok<AGROW>
            end
            R=pack(records,selected,calibration,scores,manifest,cfg);
            mri_tv_atomic_save(fullfile(cfg.outDir,'progress.mat'),'R',R);
            mri_tv_export(records,selected,calibration,manifest,cfg,false);
        end
    end
    if cfg.compareBackends
        mi=find(strcmp(cfg.methodKinds,'MRI'),1); m=selected{im,mi};
        opF=mri_tv_operators(size(xtrue),'fft_reflect');
        for r=1:nRuns
            seed=cfg.seedBase+r; b=save_observation(clean,seed,cfg,im,false);
            rng(75000+1000*im+r,'twister'); be={'spatial','fft_reflect'}; be=be(randperm(2));
            for bb=be
                if strcmp(bb{1},'spatial'), ob=opData; else, ob=opF; end
                taskid=sprintf('back_i%02d_%s_s%d',im,bb{1},seed);
                rec=task(taskid,im,seed,b,xtrue,lambdaTV,ob,m,cfg,'backend',false);
                records{end+1}=rec; %#ok<AGROW>
            end
        end
    end
end
R=pack(records,selected,calibration,scores,manifest,cfg); R.complete=true;
if cfg.tuningOnly
    R.allCalibrationSuccessful=all_selected_success(R);
    R.created=manifest.created;
    mri_tv_atomic_save(fullfile(cfg.outDir,'tuning.mat'),'T',R);
    write_scores(R);
    if ~R.allCalibrationSuccessful
        warning('MRI:TuningIncomplete', ['At least one selected method failed calibration. ' ...
            'Publication is blocked by default; inspect calibration_scores.csv.']);
    end
else
    mri_tv_atomic_save(fullfile(cfg.outDir,'benchmark.mat'),'R',R);
end
mri_tv_export(records,selected,calibration,manifest,cfg,~cfg.tuningOnly);
fprintf('\nStage %s completed: %s\n',cfg.profile,cfg.outDir);
end

function [chosen,cal,scores]=calibrate_image(im,truth,clean,lam,op,cfg)
% Randomize candidate/method order within each calibration realization.
kinds=cfg.methodKinds; nMethods=numel(kinds); lists=cell(1,nMethods);
chosen=cell(1,nMethods); scores=cell(1,nMethods); jobs=zeros(0,2); cal={};
for j=1:nMethods
    lists{j}=mri_tv_parameters(kinds{j},true,cfg.gridName);
    jobs=[jobs; [j*ones(numel(lists{j}),1),(1:numel(lists{j}))']]; %#ok<AGROW>
end
cc=cfg; cc.maxIter=cfg.calibrationMaxIter;
for r=1:numel(cfg.calibrationSeeds)
    seed=cfg.calibrationSeeds(r); b=save_observation(clean,seed,cfg,im,true);
    rng(43000+1000*im+r,'twister'); order=randperm(size(jobs,1));
    for h=order
        j=jobs(h,1); k=jobs(h,2); m=lists{j}(k);
        taskid=sprintf('cal_i%02d_%s_c%02d_s%d',im,m.kind,m.id,seed);
        cal{end+1}=task(taskid,im,seed,b,truth,lam,op,m,cc,'calibration',false); %#ok<AGROW>
    end
end
for j=1:nMethods
    cand=lists{j}; sc=zeros(numel(cand),4);
    for k=1:numel(cand)
        rr=cal(cellfun(@(z)strcmp(z.method.kind,kinds{j}) && z.method.id==cand(k).id,cal));
        assert(numel(rr)==numel(cfg.calibrationSeeds),'MRI:Calibration','Missing calibration tasks.');
        failed=cellfun(@(z)~z.success,rr); times=cellfun(@(z)z.elapsed,rr);
        gaps=cellfun(@(z)z.final.relativeGap,rr); gaps(~isfinite(gaps))=Inf;
        timeScore=Inf; if ~any(failed), timeScore=median(times); end
        sc(k,:)=[sum(failed),timeScore,median(gaps),cand(k).id];
    end
    [~,ix]=sortrows(sc,[1 2 3 4]); chosen{j}=cand(ix(1)); scores{j}=sc;
    fprintf('Selected image %d %s: candidate %d; %d/%d calibration successes.\n', ...
        im,kinds{j},chosen{j}.id,numel(cfg.calibrationSeeds)-sc(ix(1),1),numel(cfg.calibrationSeeds));
end
end

function [selected,cal,scores,info]=read_frozen(cfg,truth,lam,code)
selected={}; cal={}; scores={}; info=struct();
if isempty(cfg.tuningFile), return; end
assert(cfg.tune && ~cfg.tuningOnly,'MRI:TuningFile','Invalid frozen-tuning mode.');
if ~exist(cfg.tuningFile,'file'), error('MRI:MissingTuning','Tuning file not found: %s',cfg.tuningFile); end
z=load(cfg.tuningFile,'T'); assert(isfield(z,'T'),'MRI:TuningFile','Expected variable T.'); T=z.T;
assert(isfield(T,'complete') && T.complete && T.config.tuningOnly, ...
    'MRI:TuningFile','Incomplete/non-tuning input.');
assert(isequaln(T.manifest.truth,truth) && T.manifest.lambdaTV==lam, ...
    'MRI:TuningMismatch','Tuning image/weight differs. Run tuning for the actual image(s).');
assert(isequaln(T.manifest.commonKey,mri_tv_common_key(cfg)), ...
    'MRI:TuningMismatch','Grid, tolerances, caps, model, or solver settings differ from tuning.');
assert(isequaln(T.manifest.codeSnapshot,code), ...
    'MRI:TuningMismatch','Core source differs from frozen tuning. Use a new output root.');
assert(strcmp(T.manifest.version,version) && strcmp(T.manifest.computer,computer), ...
    'MRI:TuningMismatch','Tuning runtime differs from test runtime.');
if cfg.requireSuccessfulTuning
    assert(all_selected_success(T),'MRI:CalibrationFailed', ...
        'Selected candidates did not pass all calibration runs. Inspect tuning/calibration_scores.csv.');
end
selected=T.selected; cal=T.calibration; scores=T.scores;
info.path=cfg.tuningFile; info.created=T.created; info.allCalibrationSuccessful=all_selected_success(T);
end

function yes=all_selected_success(R)
yes=true;
for im=1:size(R.selected,1)
    for j=1:size(R.selected,2)
        m=R.selected{im,j}; if isempty(m), yes=false; return; end
        rr=R.calibration(cellfun(@(z)z.image==im && strcmp(z.method.kind,m.kind) && ...
            isequaln(z.method,m),R.calibration));
        if numel(rr)~=numel(R.config.calibrationSeeds) || any(cellfun(@(z)~z.success,rr))
            yes=false; return;
        end
    end
end
end

function b=save_observation(clean,seed,cfg,im,isCalibration)
rng(seed,'twister'); b=clean+cfg.noiseStd*randn(size(clean));
if cfg.clipObservation, b=min(255,max(0,b)); end
if isCalibration, folder='data_calibration'; else, folder='data'; end
path=fullfile(cfg.outDir,folder,sprintf('image%02d_seed%d.mat',im,seed));
if exist(path,'file')
    q=load(path,'observation');
    assert(isfield(q,'observation') && isequaln(q.observation.b,b), ...
        'MRI:DataMismatch','Saved observation mismatch: %s',path);
else
    observation=struct('b',b,'seed',seed,'noiseStd',cfg.noiseStd,'generator','twister', ...
        'clipped',cfg.clipObservation);
    mri_tv_atomic_save(path,'observation',observation);
end
end

function rec=task(key,im,seed,b,truth,lam,op,m,cfg,panel,keepState)
path=fullfile(cfg.outDir,'tasks',[key '.mat']);
if exist(path,'file') && cfg.resume
    z=load(path,'rec'); rec=z.rec;
    assert(isfield(rec,'complete') && rec.complete && isequaln(rec.method,m) && ...
        strcmp(rec.backend,op.backend) && rec.seed==seed && strcmp(rec.panel,panel), ...
        'MRI:CachedTask','Cached task does not match request: %s',key);
    fprintf('[cached] %s\n',key);
    restore_image(rec,cfg);
    return;
end
mri_tv_check_parameters(m);
if cfg.warmupIterations>0
    wc=cfg; wc.maxIter=cfg.warmupIterations; wc.checkEvery=max(1,cfg.warmupIterations);
    mri_tv_solve(b,lam,op,m,wc);
end
fprintf('%s | alpha %.4g rho %.4g tau %.4g sigma %.4g\n',key,m.alpha,m.rho,m.tau,m.sigma);
[x,u,rec]=mri_tv_solve(b,lam,op,m,cfg);
xf=min(255,max(0,x)); mse=mean((xf(:)-truth(:)).^2);
if strcmp(rec.status,'numerical_failure'), rec.psnr=NaN;
elseif mse==0, rec.psnr=Inf; else, rec.psnr=10*log10(255^2/mse); end
rec.image=im; rec.seed=seed; rec.panel=panel; rec.backend=op.backend;
rec.pixels=numel(b); rec.lambdaTV=lam; rec.key=key;
rec.parameterCheck=mri_tv_check_parameters(m); rec.imagePath='';
rec.xRaw=[]; rec.uRaw=[];
if keepState
    rec.xRaw=x; rec.uRaw=u;
    rec.imagePath=sprintf('i%02d_%s.png',im,m.kind);
end
rec.complete=true;
mri_tv_atomic_save(path,'rec',rec); restore_image(rec,cfg);
fprintf('  %s | %d iterations | %.3f s | KKT %.2e | gap %.2e | feas %.2e\n', ...
    rec.status,rec.iter,rec.elapsed,rec.final.kkt,rec.final.relativeGap,rec.final.feas);
end

function restore_image(rec,cfg)
if cfg.saveImages && ~isempty(rec.imagePath) && ~isempty(rec.xRaw) && all(isfinite(rec.xRaw(:)))
    imwrite(uint8(min(255,max(0,rec.xRaw))),fullfile(cfg.outDir,'figures',rec.imagePath));
end
end

function grid=sensitivity_grid(cfg)
grid={}; id=0;
for rho=cfg.sensitivityRho
    for a=cfg.sensitivityAlpha
        for kk={'MRI','RI'}
            m=struct('kind',kk{1},'alpha',a,'rho',rho,'tau',cfg.sensitivityTau, ...
                'sigma',cfg.sensitivitySigma,'id',0);
            try
                mri_tv_check_parameters(m);
            catch ex
                if strcmp(ex.identifier,'MRI:Inadmissible'), continue; end
                rethrow(ex);
            end
            id=id+1; m.id=id; grid{end+1}=m; %#ok<AGROW>
        end
    end
end
end

function R=pack(records,selected,calibration,scores,manifest,cfg)
R.records=records; R.selected=selected; R.calibration=calibration;
R.scores=scores; R.manifest=manifest; R.config=cfg; R.complete=false;
end

function write_grid(cfg)
fid=fopen(fullfile(cfg.outDir,'candidate_parameters.csv'),'w');
assert(fid>=0,'MRI:Output','Cannot create candidate_parameters.csv.'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Method,Candidate,Alpha,Rho,Tau,Sigma,OperatorBudget\n');
for j=1:numel(cfg.methodKinds)
    ms=mri_tv_parameters(cfg.methodKinds{j},true,cfg.gridName);
    for m=ms
        a=mri_tv_check_parameters(m);
        fprintf(fid,'%s,%d,%.17g,%.17g,%.17g,%.17g,%.17g\n', ...
            m.kind,m.id,m.alpha,m.rho,m.tau,m.sigma,a.operatorBudget);
    end
end
end

function write_scores(R)
fid=fopen(fullfile(R.config.outDir,'calibration_scores.csv'),'w');
assert(fid>=0,'MRI:Output','Cannot create calibration_scores.csv.'); c=onCleanup(@()fclose(fid)); %#ok<NASGU>
fprintf(fid,'Image,Method,Candidate,Failed,MedianTimeScore,MedianRelativeGap,Selected\n');
for im=1:size(R.scores,1)
    for j=1:size(R.scores,2)
        s=R.scores{im,j}; m=R.selected{im,j};
        for k=1:size(s,1)
            fprintf(fid,'%d,%s,%d,%d,%.17g,%.17g,%d\n',im,m.kind,s(k,4),s(k,1),s(k,2),s(k,3),s(k,4)==m.id);
        end
    end
end
end

function release_lock(path)
if exist(path,'dir'), rmdir(path); end
end

function validate_cfg(c,nRuns,lam)
validateattributes(lam,{'numeric'},{'scalar','real','finite','positive'});
validateattributes(nRuns,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
assert(nRuns>0 || c.tuningOnly,'MRI:Runs','Zero test runs are only allowed for tuning.');
for nm={'maxIter','calibrationMaxIter','checkEvery'}
    validateattributes(c.(nm{1}),{'numeric'},{'scalar','real','finite','integer','positive'});
end
for nm={'tolKKT','tolGap','tolFeas','diagnosticTau','diagnosticSigma'}
    validateattributes(c.(nm{1}),{'numeric'},{'scalar','real','finite','positive'});
end
validateattributes(c.noiseStd,{'numeric'},{'scalar','real','finite','nonnegative'});
validateattributes(c.warmupIterations,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
validateattributes(c.calibrationSeeds,{'numeric'},{'vector','integer','nonnegative','finite','nonempty'});
validateattributes(c.seedBase,{'numeric'},{'scalar','integer','nonnegative','finite'});
assert(numel(unique(c.calibrationSeeds))==numel(c.calibrationSeeds), ...
    'MRI:Seeds','Duplicate calibration seeds.');
assert(~any(ismember(c.calibrationSeeds,c.seedBase+(1:nRuns))), ...
    'MRI:Seeds','Calibration and test seeds overlap.');
assert(numel(unique(c.methodKinds))==numel(c.methodKinds) && ...
    all(ismember(c.methodKinds,{'MRI','RI','CV','CP'})) && ~isempty(c.methodKinds), ...
    'MRI:Methods','Method list is empty, duplicated, or invalid.');
if c.compareBackends
    assert(any(strcmp(c.methodKinds,'MRI')),'MRI:Methods','Backend comparison requires MRI.');
end
assert(~c.tuningOnly || c.tune,'MRI:Mode','Tuning-only requires tuning enabled.');
end
