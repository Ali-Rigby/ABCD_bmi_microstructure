function vols = process_spline_derivatives(vols, tbl_bf, xvec, beta_hat, idx_bf, idx_intercept, mask, coeffCovar, slice_indices, modality, designmat_label, colnames_model, limits, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, xname, vol_display_mode, gray_limits)
%PROCESS_SPLINE_DERIVATIVES Append spline-derivative maps to a showVol volume cell array.
%
% Uses the 1st output of FEMA_convert_splinevols (vols_dbf: 4D derivative volumes).
% Note: the 4th output is volmat_dbf (2D voxel matrix), not a volume — do not use it here.
% Each slice gets grayscale + colormap entries via convert_and_append_vols
% vol_display_mode: 'both' | 'color' | 'gray' (default 'both').

    if nargin < 21 || isempty(vol_display_mode)
        vol_display_mode = 'both';
    end
    if nargin < 22 || isempty(gray_limits)
        gray_limits = [];
    end

    if isempty(idx_bf)
        return;
    end

    vols_deriv = FEMA_convert_splinevols( ...
        tbl_bf, xvec, beta_hat, idx_bf, idx_intercept, mask, coeffCovar);

    n_slices = size(vols_deriv, 4);
    if n_slices <= 1 && ndims(vols_deriv) < 4
        n_slices = 1;
    end

    slice_indices = slice_indices(:);
    for k = 1:numel(slice_indices)
        slice = slice_indices(k);
        if slice < 1 || slice > n_slices
            continue;
        end
        statname = sprintf('beta_hat_derivative (%s = %.1f), %s', xname, xvec(slice), designmat_label);
        vols = convert_and_append_vols(vols, vols_deriv(:, :, :, slice), modality, statname, ...
            colnames_model(idx_bf), 1, limits, showPval, cmap, bg_vol, CSFprob, ...
            interpMethod, atlasVersion, vol_display_mode, gray_limits);
    end
end
