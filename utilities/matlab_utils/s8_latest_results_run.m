function run_name = s8_latest_results_run(results_root)
%S8_LATEST_RESULTS_RUN Newest existing results_YYYY_MM_DD under results_root.

run_name = '';
if ~exist(results_root, 'dir')
    return;
end

d = dir(fullfile(results_root, 'results_*'));
d = d([d(:).isdir]);
if isempty(d)
    return;
end

names = {d(:).name};
[~, ord] = sort(names, 'descend');
run_name = names{ord(1)};

end
