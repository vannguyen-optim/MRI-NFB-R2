function snapshot = mri_tv_source_snapshot()
% Record exact .m source in the core folder, excluding reporting/entry scripts.
folder=fileparts(mfilename('fullpath')); ff=dir(fullfile(folder,'*.m'));
[~,idx]=sort({ff.name}); ff=ff(idx); snapshot=cell(numel(ff),2);
for j=1:numel(ff)
    snapshot{j,1}=ff(j).name;
    snapshot{j,2}=fileread(fullfile(folder,ff(j).name));
end
end
