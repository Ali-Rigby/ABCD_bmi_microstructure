function T = bmi_rank_roi_wald(vol_W, brain_mask, cfg, min_voxels, vol_FDR)
%BMI_RANK_ROI_WALD Rank ROIs by mean uncorrected Wald statistic (descending).
%
% Uses s7 Wald volumes (vol_W_s etc.) and atlas definitions from bmi_roi_wald_config.
% All voxels in each ROI (within the FEMA brain mask) contribute to the mean.
% Optional vol_FDR (vol_pW_*_FDR): adds n_fdr_sig = count of FDR-significant voxels in ROI.

if nargin < 4 || isempty(min_voxels)
    min_voxels = 1;
end
if nargin < 5
    vol_FDR = [];
end

rows = {};
row_i = 0;

aseg = bmi_load_atlas_context(cfg.dirname_atlas, 'aseg');
fiber = bmi_load_atlas_context(cfg.dirname_atlas, 'fiber');
pauli = bmi_load_atlas_context(cfg.dirname_atlas, 'pauli_subcortnuclei_abcd3_cor10');

% --- Bilateral aseg ---
for ii = 1:numel(cfg.bilateral_aseg)
    base = cfg.bilateral_aseg{ii};
    left_mask = bmi_build_roi_mask(aseg.atlas_type, aseg.roinames, aseg.volmat_prob, ...
        {['Left-' base]}, cfg.prob_thresh);
    right_mask = bmi_build_roi_mask(aseg.atlas_type, aseg.roinames, aseg.volmat_prob, ...
        {['Right-' base]}, cfg.prob_thresh);
    if ~any(left_mask(:)) && ~any(right_mask(:))
        continue;
    end
    vol_roi = left_mask | right_mask;
    stats = bmi_roi_mean_wald_from_volume(vol_W, brain_mask, vol_roi, min_voxels, vol_FDR);
    row_i = row_i + 1;
    rows(row_i, :) = {'bilateral', ['Bilateral-' base], stats.mean_wald, stats.median_wald, ...
        stats.n_voxels, stats.n_fdr_sig}; %#ok<AGROW>
end

% --- Fiber (single) ---
for ii = 1:numel(cfg.fiber_single)
    name = cfg.fiber_single{ii};
    vol_roi = bmi_build_roi_mask(fiber.atlas_type, fiber.roinames, fiber.volmat_prob, {name}, cfg.prob_thresh);
    stats = bmi_roi_mean_wald_from_volume(vol_W, brain_mask, vol_roi, min_voxels, vol_FDR);
    row_i = row_i + 1;
    rows(row_i, :) = {'fiber', name, stats.mean_wald, stats.median_wald, ...
        stats.n_voxels, stats.n_fdr_sig}; %#ok<AGROW>
end

% --- Fiber (bilateral L/R) ---
for ii = 1:numel(cfg.fiber_bilateral)
    base = cfg.fiber_bilateral{ii};
    left_mask = bmi_build_roi_mask(fiber.atlas_type, fiber.roinames, fiber.volmat_prob, ...
        {['L_' base]}, cfg.prob_thresh);
    right_mask = bmi_build_roi_mask(fiber.atlas_type, fiber.roinames, fiber.volmat_prob, ...
        {['R_' base]}, cfg.prob_thresh);
    vol_roi = left_mask | right_mask;
    stats = bmi_roi_mean_wald_from_volume(vol_W, brain_mask, vol_roi, min_voxels, vol_FDR);
    row_i = row_i + 1;
    rows(row_i, :) = {'fiber', ['Bilateral-' base], stats.mean_wald, stats.median_wald, ...
        stats.n_voxels, stats.n_fdr_sig}; %#ok<AGROW>
end

% --- Pauli ---
for ii = 1:numel(cfg.pauli)
    name = cfg.pauli{ii};
    vol_roi = bmi_build_roi_mask(pauli.atlas_type, pauli.roinames, pauli.volmat_prob, {name}, cfg.prob_thresh);
    stats = bmi_roi_mean_wald_from_volume(vol_W, brain_mask, vol_roi, min_voxels, vol_FDR);
    row_i = row_i + 1;
    rows(row_i, :) = {'pauli', name, stats.mean_wald, stats.median_wald, ...
        stats.n_voxels, stats.n_fdr_sig}; %#ok<AGROW>
end

if row_i == 0
    T = table();
    return;
end

T = cell2table(rows, 'VariableNames', ...
    {'roi_group', 'roi_name', 'mean_wald', 'median_wald', 'n_voxels', 'n_fdr_sig'});
T = sortrows(T, 'mean_wald', 'descend');
T.rank = (1:height(T))';

% Reorder columns
T = T(:, {'rank', 'roi_group', 'roi_name', 'mean_wald', 'median_wald', 'n_voxels', 'n_fdr_sig'});

end
