function stats = bmi_roi_mean_wald_from_volume(vol_W, brain_mask, vol_roi, min_voxels, vol_FDR)
%BMI_ROI_MEAN_WALD_FROM_VOLUME Mean Wald statistic within an ROI mask.
%
% Averages uncorrected Wald stats over all in-mask voxels in the ROI (not
% FDR- or Bonferroni-filtered). ROI masks are subsampled (2x2x2) to match
% atlas space; vol_W and brain_mask are subsampled when sizes differ.
%
% Optional vol_FDR (same space as vol_W): s7 FDR indicator volume where
% 0 = FDR-significant and 1 = not. When provided, stats.n_fdr_sig counts
% ROI ∩ brain voxels that are FDR-significant.
%
% Inputs:
%   vol_W       - 3-D Wald statistic volume (s7 vol_W_s / vol_W_l / vol_W_ls)
%   brain_mask  - FEMA analysis mask (same space as vol_W before subsampling)
%   vol_roi     - logical 3-D ROI mask (subsampled atlas space)
%   min_voxels  - minimum voxels required (default 1)
%   vol_FDR     - optional FDR volume (vol_pW_*_FDR); 0 = significant
%
% Output struct fields: mean_wald, n_voxels, median_wald, n_fdr_sig

if nargin < 4 || isempty(min_voxels)
    min_voxels = 1;
end
if nargin < 5
    vol_FDR = [];
end

vol_roi = logical(vol_roi);

if ~isequal(size(vol_roi), size(vol_W))
    vol_W = subsample_volume(vol_W);
    brain_mask = subsample_volume(brain_mask >= 0.5);
    if ~isempty(vol_FDR)
        vol_FDR = subsample_volume(vol_FDR);
    end
else
    brain_mask = brain_mask >= 0.5;
end

if ~isequal(size(vol_roi), size(vol_W))
    error('bmi_roi_mean_wald_from_volume: ROI mask size [%s] != Wald volume size [%s] after subsampling.', ...
        num2str(size(vol_roi)), num2str(size(vol_W)));
end

combined = vol_roi & brain_mask;
vals = vol_W(combined);
vals = vals(isfinite(vals));

stats = struct('mean_wald', NaN, 'median_wald', NaN, 'n_voxels', 0, 'n_fdr_sig', 0);
if numel(vals) < min_voxels
    return;
end

stats.mean_wald = mean(vals);
stats.median_wald = median(vals);
stats.n_voxels = numel(vals);

if ~isempty(vol_FDR)
    if ~isequal(size(vol_FDR), size(vol_W))
        error('bmi_roi_mean_wald_from_volume: FDR volume size mismatch after subsampling.');
    end
    % s7 convention: 0 = FDR significant, 1 = not
    stats.n_fdr_sig = nnz(combined & (vol_FDR == 0));
end

end
