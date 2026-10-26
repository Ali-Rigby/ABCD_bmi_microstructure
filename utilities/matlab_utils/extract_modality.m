function modality = extract_modality(fname_results)
%EXTRACTMODALITY Imaging modality tag (FEMA fstem_imaging) from results path.
%
% Examples:
%   .../FEMA_wrapper_output_voxel_rni.mat           -> 'rni'
%   .../FEMA_wrapper_output_voxel_rni_estimates.mat -> 'rni'

[~, name, ~] = fileparts(fname_results);
[~, fstem_imaging] = fema_parse_results_stem(name);
modality = fstem_imaging;
end
