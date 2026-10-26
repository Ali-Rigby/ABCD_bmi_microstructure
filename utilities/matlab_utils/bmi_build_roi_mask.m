function vol_roi = bmi_build_roi_mask(atlas_type, roinames, volmat_prob, roi_names, prob_thresh)
%BMI_BUILD_ROI_MASK Binary ROI mask from atlas probability volume(s).
%
% roi_names - cellstr of label fragments to match (via contains), combined with OR.

if nargin < 5 || isempty(prob_thresh)
    prob_thresh = 0.7;
end

idx_roi = [];
for ii = 1:numel(roi_names)
    hit = find(contains(roinames, roi_names{ii}));
    if isempty(hit)
        warning('bmi_build_roi_mask: ''%s'' not found in atlas %s', roi_names{ii}, atlas_type);
    else
        idx_roi = [idx_roi, hit(:)']; %#ok<AGROW>
    end
end

idx_roi = unique(idx_roi, 'stable');
if isempty(idx_roi)
    vol_roi = false(size(volmat_prob, 1), size(volmat_prob, 2), size(volmat_prob, 3));
    return;
end

switch atlas_type
    case {'aseg', 'fiber'}
        % volmat_prob: nroi x dim1 x dim2 x dim3 (MATLAB load) or permuted
        if size(volmat_prob, 1) == numel(roinames)
            prob_stack = squeeze(max(volmat_prob(idx_roi, :, :, :), [], 1));
        else
            prob_stack = squeeze(max(volmat_prob(:, :, :, idx_roi), [], 4));
        end
    case {'thalamus_abcd3_cor10', 'pauli_subcortnuclei_abcd3_cor10'}
        % vol_labels_reg: dim1 x dim2 x dim3 x nroi
        prob_stack = max(volmat_prob(:, :, :, idx_roi), [], 4);
    otherwise
        if ndims(volmat_prob) == 4 && size(volmat_prob, 4) == numel(roinames)
            prob_stack = max(volmat_prob(:, :, :, idx_roi), [], 4);
        elseif size(volmat_prob, 1) == numel(roinames)
            prob_stack = squeeze(max(volmat_prob(idx_roi, :, :, :), [], 1));
        else
            error('bmi_build_roi_mask: unrecognized volmat_prob shape for %s', atlas_type);
        end
end

if ndims(prob_stack) == 4
    prob_stack = squeeze(prob_stack);
end

vol_roi = prob_stack > prob_thresh;
vol_roi = subsample_volume(vol_roi);
vol_roi = logical(vol_roi);

end
