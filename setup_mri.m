function root = setup_mri()
%SETUP_MRI Add only this distribution to the beginning of the MATLAB path.
root=fileparts(mfilename('fullpath'));
addpath(root,'-begin');
addpath(fullfile(root,'matlab'),'-begin');
addpath(fullfile(root,'reporting'),'-begin');
addpath(fullfile(root,'tests'),'-begin');
rehash;
folders={root,fullfile(root,'matlab'),fullfile(root,'reporting'),fullfile(root,'tests')};
for j=1:numel(folders)
    entries=dir(fullfile(folders{j},'*.m'));
    for k=1:numel(entries)
        [~,name]=fileparts(entries(k).name);
        actual=which(name); expected=fullfile(folders{j},entries(k).name);
        assert(strcmpi(actual,expected),'MRI:PathShadow', ...
            'Another file shadows %s. Change Current Folder to this package root: %s',name,root);
    end
end
end
