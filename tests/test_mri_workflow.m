function checks = test_mri_workflow()
%TEST_MRI_WORKFLOW Tiny constant-image integration test; NOT research data.
% Exercises tuning, frozen import, per-task resume, export, and refusal guards.
root=setup_mri(); folder=tempname; mkdir(folder);
guard=onCleanup(@()cleanup_dir(folder)); %#ok<NASGU>
image=100*ones(4,5); cfg=mri_tv_config('tuning');
cfg.outDir=fullfile(folder,'tuning'); cfg.runTests=false;
cfg.noiseStd=0;
cfg.maxIter=3; cfg.calibrationMaxIter=3; cfg.warmupIterations=0;
cfg.makePlots=false; cfg.saveImages=false;
text=evalc('T=run_mri_tv(image,3,0,cfg);');
assert(T.complete && T.allCalibrationSuccessful);
nCandidates=numel(mri_tv_parameters('MRI',true,cfg.gridName));
assert(numel(T.calibration)==4*nCandidates*numel(cfg.calibrationSeeds));
assert(numel(T.selected)==4);
text=[text evalc('U=run_mri_tv(image,3,0,cfg);')]; %#ok<AGROW>
assert(isequaln(T.selected,U.selected));
assert(isequaln(T.calibration,U.calibration));
c=cfg; c.profile='publication'; c.tuningOnly=false;
c.outDir=fullfile(folder,'test'); c.tuningFile=fullfile(folder,'tuning','tuning.mat');
text=[text evalc('R=run_mri_tv(image,3,2,c);')]; %#ok<AGROW>
assert(numel(R.records)==8 && all(cellfun(@(r)r.success,R.records)));
text=[text evalc('R2=run_mri_tv(image,3,2,c);')]; %#ok<AGROW>
assert(isequaln(R.records,R2.records));
rep=export_mri_tv_paper(c.outDir,false); %#ok<NASGU>
assert(exist(fullfile(c.outDir,'summary_benchmark.csv'),'file')==2);
t=fileread(fullfile(c.outDir,'tables','paper_parameters_i01.tex'));
assert(~isempty(strfind(t,'\begin{tabular}'))); %#ok<STREMP>
assert(~any(uint8(t)<32 & ~ismember(uint8(t),[9 10 13])));
changed=c; changed.tolKKT=2*c.tolKKT;
expect_error(@()run_mri_tv(image,3,2,changed),'MRI:TuningMismatch');
expect_error(@()run_mri_tv(image+1,3,2,c),'MRI:TuningMismatch');
assert(exist(fullfile(c.outDir,'.running_lock'),'dir')==0);
bad=T;
m=bad.selected{1,1}; ix=find(cellfun(@(r)strcmp(r.method.kind,m.kind) && isequaln(r.method,m),bad.calibration),1);
bad.calibration{ix}.success=false;
badfile=fullfile(folder,'failed_tuning.mat'); mri_tv_atomic_save(badfile,'T',bad);
changed=c; changed.tuningFile=badfile;
expect_error(@()run_mri_tv(image,3,2,changed),'MRI:CalibrationFailed');
checks=struct('passed',true,'calibrationTasks',numel(T.calibration),'testTasks',numel(R.records), ...
    'resumeIdentical',true,'tuningMismatchGuard',true,'failedCalibrationGuard',true,'tableExport',true);
fid=fopen(fullfile(root,'validation','workflow_test_log.txt'),'w');
assert(fid>=0); fprintf(fid,'%s',text); fclose(fid);
end
function expect_error(f,id)
received='';
try, f(); catch ex, received=ex.identifier; end
assert(strcmp(received,id),'MRI:GuardTest','Expected %s; received %s.',id,received);
end
function cleanup_dir(folder)
if exist(folder,'dir'), rmdir(folder,'s'); end
end
