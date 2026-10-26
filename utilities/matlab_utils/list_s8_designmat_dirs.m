function names = list_s8_designmat_dirs(dirname_out)
%LIST_S8_DESIGNMAT_DIRS Designmat result subfolders under a results_* directory.
d = dir(fullfile(dirname_out, 'designmat*'));
d = d([d.isdir]);
names = {d.name};
keep = false(size(names));
for i = 1:numel(names)
    % Sex-split (1f/1m) or combined-sex (designmat1_...)
    keep(i) = ~isempty(regexp(names{i}, '^designmat\d+[fm]?_', 'once'));
end
names = sort(names(keep));
end
