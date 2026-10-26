function fpath = fema_results_path(dirname_out, datatype, fstem_imaging)
%FEMA_RESULTS_PATH Resolve FEMA 4.x / legacy estimates .mat for one run.
%
% Search order:
%   1. FEMA_wrapper_output_<datatype>_<fstem>_estimates.mat  (s6 outPrefix)
%   2. FEMA_wrapper_output_<datatype>_<fstem>.mat           (legacy single file)
%   3. FEMA_estimates.mat                                   (bare 4.x, one run per folder)

pref = sprintf('FEMA_wrapper_output_%s_%s', datatype, fstem_imaging);
candidates = {
    fullfile(dirname_out, [pref, '_estimates.mat'])
    fullfile(dirname_out, [pref, '.mat'])
    fullfile(dirname_out, 'FEMA_estimates.mat')
    };

for i = 1:numel(candidates)
    if isfile(candidates{i})
        fpath = candidates{i};
        return;
    end
end

error('fema_results_path: no estimates file under %s for %s / %s.', ...
    dirname_out, datatype, fstem_imaging);
end
