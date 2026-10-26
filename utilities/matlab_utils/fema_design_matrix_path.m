function fpath = fema_design_matrix_path(dirname_out, datatype, fstem_imaging)
%FEMA_DESIGN_MATRIX_PATH Resolve FEMA 4.x design-matrix .mat for one run.
%
% Search order:
%   1. FEMA_wrapper_output_<datatype>_<fstem>_designMatrix.mat  (s6 outPrefix)
%   2. FEMA_designMatrix.mat
%   3. FEMA_wrapper_output_<datatype>_<fstem>.mat                 (legacy single file)

pref = sprintf('FEMA_wrapper_output_%s_%s', datatype, fstem_imaging);
candidates = {
    fullfile(dirname_out, [pref, '_designMatrix.mat'])
    fullfile(dirname_out, 'FEMA_designMatrix.mat')
    fullfile(dirname_out, [pref, '.mat'])
    };

for i = 1:numel(candidates)
    if isfile(candidates{i})
        fpath = candidates{i};
        return;
    end
end

error('fema_design_matrix_path: no design-matrix file under %s for %s / %s.', ...
    dirname_out, datatype, fstem_imaging);
end
