function vols = process_spline_slices(vols, tbl_bf, xvec, beta_hat, idx_bf, idx_intercept, mask, coeffCovar, slice_indices, modality, designmat_label, colnames_model, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, xname)

% process_spline_slices  Process spline basis function slices and append voxelwise volumes.
%
% Usage:
%   vols = process_spline_slices(vols, tbl_bf, xvec, beta_hat, idx_bf, ...
%                                 idx_intercept, mask, coeffCovar, slice_indices, ...
%                                 modality, designmat_label, colnames_model, limits, ...
%                                 showPval, cmap, bg_vol, CSFprob, interpMethod, ...
%                                 atlasVersion, xname)
%
% Inputs:
%   vols             - Existing cell array of volumes (can be empty) to append results to.
%   tbl_bf           - Table of basis function values used in the model.
%   xvec             - Vector of the continuous variable modeled with splines.
%   beta_hat         - FEMA beta estimates.
%   idx_bf           - Indices of basis function terms in the model.
%   idx_intercept    - Index of intercept term in the model.
%   mask             - Voxel mask for analysis.
%   coeffCovar       - Covariance matrix of model coefficients.
%   slice_indices    - Vector of indices into xvec indicating which slices to extract.
%   modality         - String label for imaging modality (e.g., 'rni').
%   designmat_label  - String label for the design matrix (for figure labeling).
%   colnames_model   - Cell array of column names from the design matrix.
%   limits           - Color scale limits for visualization (e.g., [-0.01 0.01]).
%   showPval         - Logical flag to overlay p-values instead of beta estimates.
%   cmap             - Colormap to use for visualization.
%   bg_vol           - Background anatomical volume for overlay.
%   CSFprob          - CSF probability threshold for masking.
%   interpMethod     - Method for upsampling ('linear' or 'nearest').
%   atlasVersion     - String specifying atlas version (e.g., 'ABCD3').
%   xname            - Name of the predictor variable (for labeling).
%
% Output:
%   vols             - Updated cell array of volumes with appended slices for each spline value.
%
% Description:
%   This function generates voxelwise maps for selected spline slices of a continuous predictor
%   variable. For each specified slice index, it extracts the corresponding beta map, converts
%   it to a 1mm resolution volume using the provided settings, and appends it to the input `vols`.
%
%   Each appended volume will have a name in the format:
%       'beta_hat (xname = value), designmat_label'
%
% Example:
%   vols = process_spline_slices([], tbl_bf_bmi, bmivec, beta_hat, Sbmi_idx, ...
%          intercept_idx, mask, coeffCovar, [2, 11, 21], 'rni', 'designmat1', ...
%          colnames_model, [-0.01 0.01], false, redblackblue, [], [], [], 'ABCD3', 'bmi');
%
% See also: FEMA_convert_splinevols, convert_and_append_vols

    if isempty(idx_bf)
        return;
    end
    [vols_tmp] = FEMA_convert_splinevols(tbl_bf, xvec, beta_hat, idx_bf, idx_intercept, mask, coeffCovar);
    for slice = slice_indices
        statname = sprintf('beta_hat (%s = %.1f), %s', xname, xvec(slice), designmat_label);
        vols = convert_and_append_vols(vols, vols_tmp(:,:,:,slice), modality, statname, colnames_model(idx_bf), 1, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion);
    end
end