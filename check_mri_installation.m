function checks = check_mri_installation(runConvergence)
%CHECK_MRI_INSTALLATION Static diagnostics plus executable mathematical tests.
% No optimization/image-processing/parallel toolboxes are needed.
if nargin<1, runConvergence=true; end
root=setup_mri(); reports={};
if exist('checkcode','file')
    folders={root,fullfile(root,'matlab'),fullfile(root,'reporting'),fullfile(root,'tests')};
    for f=1:numel(folders)
        ff=dir(fullfile(folders{f},'*.m'));
        for j=1:numel(ff)
            file=fullfile(folders{f},ff(j).name);
            info=checkcode(file,'-id');
            reports{end+1}=struct('file',file,'issues',info); %#ok<AGROW>
            % PARSE diagnostics are fatal; other Code Analyzer advice is saved.
            if ~isempty(info) && isfield(info,'id')
                ids={info.id};
                if any(strcmp(ids,'PARSE'))
                    error('MRI:Parse','MATLAB parser error in %s.',file);
                end
            end
        end
    end
else
    warning('MRI:Analyzer','MATLAB checkcode is unavailable; executable tests still run.');
end
checks=test_mri_tv(runConvergence);
if runConvergence, checks.workflow=test_mri_workflow(); end
checks.sourceFiles=mri_tv_source_snapshot();
checks.codeAnalyzer=reports; checks.MATLABVersion=version; checks.platform=computer;
mri_tv_atomic_save(fullfile(root,'validation','matlab_checks.mat'),'checks',checks);
end
