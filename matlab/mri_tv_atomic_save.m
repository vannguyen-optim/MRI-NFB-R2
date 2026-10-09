function mri_tv_atomic_save(path,name,value)
% Write new content to a temporary MAT file before replacing the target.
parent=fileparts(path); if ~exist(parent,'dir'), mkdir(parent); end
part=[tempname(parent) '.mat']; payload=struct(); payload.(name)=value;
try
    save(part,'-struct','payload','-v7');
    [ok,msg]=movefile(part,path,'f');
    if ~ok, error('MRI:SaveFailed','%s',msg); end
catch ex
    if exist(part,'file'), delete(part); end
    rethrow(ex);
end
end
