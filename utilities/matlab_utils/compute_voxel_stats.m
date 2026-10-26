function compute_voxel_stats(beta_hat, coeffCovar, mask, colnames_model, dirname_out, ...
                           linear_term, xname, fname_results, dirname_desmat, bf_dir, bf_dir_shared)
% Compute voxelwise linear/spline/combined stats, using BF tables from bf_dir.
%
% Inputs:
%   - beta_hat, coeffCovar, mask, colnames_model : FEMA outputs
%   - dirname_out      : where to save outputs
%   - linear_term      : e.g., 'ab_g_dyn__visit_age'
%   - xname            : e.g., 'age' | 'bmi' | 'pds'
%   - fname_results    : full path to FEMA estimates .mat (see fema_results_path)
%   - dirname_desmat   : design-matrix folder name (used to pick male/female BF for BMI)
%   - bf_dir           : profile datafiles (pheno BF tables from s3)
%   - bf_dir_shared    : shared/datafiles (age + PDS BF tables from s3); defaults to bf_dir

if nargin < 11 || isempty(bf_dir_shared)
    bf_dir_shared = bf_dir;
end
%
% Saves:
%   <modality>_Wald_splines_CI_output_<xname>.mat  (in dirname_out)


    % Parse modality from result file
    modality = extract_modality(fname_results);  
    if isempty(modality)
        warning('Could not determine modality from filename: %s', fname_results);
        modality = 'unknown';
    else
        modality = char(modality);
    end
    fprintf('Processing Voxel Data: %s (Modality: %s)\n', fname_results, modality);


    % Indices for linear, spline, and intercept
    idx_linear    = find(strcmp(colnames_model, linear_term));
    idx_splines   = find(contains(colnames_model, ['bf_', xname]));
    idx_intercept = find(contains(colnames_model, 'intercept'));

    % Multiple-comparison bookkeeping
    nvox                = numel(find(mask >= 0.5));  % # voxels within mask
    correction_factor   = nvox;
    min_pval_threshold  = 1e-300;                    % avoid log10(0)

    % Load basis functions from bf_dir
    switch lower(xname)
        case 'age'
            tbl_bf = readtable(fullfile(bf_dir_shared, 'table_bf_age.txt'), 'Delimiter','\t');
            % xvec = 8.3:0.1:17.8; % (keep if you later plot) 
        case 'bmi'
            tbl_bf = readtable(fullfile(bf_dir, 'table_bf_bmi_perc.txt'), 'Delimiter','\t');
        case 'perc_median'
            tbl_bf = readtable(fullfile(bf_dir, 'table_bf_perc_median.txt'), 'Delimiter','\t');
        case 'adj_perc_median'
            tbl_bf = readtable(fullfile(bf_dir, 'table_bf_adj_perc_median.txt'), 'Delimiter','\t');
        % case 'bmi'
        %     if contains(dirname_desmat, 'females')
        %         tbl_bf = readtable(fullfile(bf_dir, 'table_bf_bmi_perc_females.txt'), 'Delimiter','\t');
        %     elseif contains(dirname_desmat, 'males')
        %         tbl_bf = readtable(fullfile(bf_dir, 'table_bf_bmi_perc_males.txt'), 'Delimiter','\t');
        %     else
        %         error('Design matrix directory must include "females" or "males" for BMI BF tables.');
        %     end
        case 'pds'
            tbl_bf = readtable(fullfile(bf_dir_shared, 'table_bf_pds.txt'), 'Delimiter','\t');
        otherwise
            error('Unknown xname "%s" for BF loading.', xname);
    end
    bfmat = table2array(tbl_bf);

    % Flags for presence in design
    has_linear  = ~isempty(idx_linear);
    has_splines = ~isempty(idx_splines);

    % ===========================
    % LINEAR EFFECT ONLY
    % ===========================
    if has_linear
        fprintf('Performing linear effect analysis for: %s (Modality: %s)\n', linear_term, modality);

        % Linear + intercept
        statmat_linear = beta_hat([idx_linear, idx_intercept], :);
        % point-estimate at weight [1; 1] (linear + intercept)
        wvec            = [1; 1];
        volmat_linear   = statmat_linear(1, :) + statmat_linear(2, :);

        % SE for linear + intercept
        temp_sevec    = pagemtimes(coeffCovar([idx_linear, idx_intercept], [idx_linear, idx_intercept], :), wvec);
        sevec_linear  = sqrt(pagemtimes(temp_sevec, permute(wvec, [2 1])));
        semat_linear  = squeeze(sevec_linear);  % ensure [1 x V]

        % Wald test (linear term)
        L = zeros(size(beta_hat, 1));
        L(idx_linear, idx_linear) = 1;
        hypValue = 0; doF = false; numObs = size(beta_hat, 2);
        [W_l, logp_l, ~, ~, dfW_l, LW_l] = FEMA_WaldTest(L, beta_hat, coeffCovar, hypValue, doF, numObs);
        pW_l = 10.^(-logp_l);
        pW_l(isinf(logp_l)) = realmin;

        % Expand to volume
        vol_W_l      = fullvol(W_l,     mask);
        vol_pW_l     = fullvol(pW_l,    mask);
        vol_logp_l   = fullvol(logp_l,  mask);

        % Bonferroni
        vol_pW_l_Bon          = min(max(vol_pW_l * correction_factor, min_pval_threshold), 1);
        vol_pW_l_Bon(mask<.5) = 1;
        vol_logp_l_Bon        = -log10(vol_pW_l_Bon);

        % FDR correction
        pvals_masked = pW_l;
        [h_FDR, crit_p] = fdr_bh(pvals_masked, 0.05, 'pdep');
        vol_pW_l_FDR = ones(size(mask));
        vol_pW_l_FDR(mask >= 0.5) = double(~h_FDR);  % 1 = not significant, 0 = significant
        % For visualization: use uncorrected -log10(p) where FDR is significant, 0 elsewhere
        vol_logp_l_FDR = vol_logp_l;  % Start with uncorrected log p-values
        vol_logp_l_FDR(vol_pW_l_FDR > 0.5) = 0;  % Set non-significant regions to 0
    end

    % ===========================
    % SPLINE EFFECT ONLY
    % ===========================
    if has_splines
        fprintf('Performing spline effect analysis for: %s (Modality: %s)\n', xname, modality);

        statmat_spline = beta_hat([idx_splines, idx_intercept], :);
        % preallocate
        volmat_spline = zeros(size(bfmat,1), size(beta_hat,2), 'like', beta_hat);
        semat_spline  = zeros(size(bfmat,1), size(beta_hat,2), 'like', beta_hat);

        for xi = 1:size(bfmat, 1)
            wvec        = bfmat(xi, :)';
            valvec_bf   = sum(statmat_spline .* [wvec; 1], 1);
            sevec_bf    = sqrt( ...
                              pagemtimes( ...
                                  pagemtimes([wvec; 1]', coeffCovar([idx_splines idx_intercept], [idx_splines idx_intercept], :)), ...
                                  [wvec; 1] ...
                              ) ...
                           );
            volmat_spline(xi, :) = valvec_bf;
            semat_spline(xi, :)  = sevec_bf;
        end

        % Wald test: all spline coeffs jointly ≠ 0
        L = zeros(size(beta_hat, 1));
        L(idx_splines, idx_splines) = eye(length(idx_splines));

        hypValue = 0; doF = false; numObs = size(beta_hat, 2);
        [W_s, logp_s, ~, ~, dfW_s, LW_s] = FEMA_WaldTest(L, beta_hat, coeffCovar, hypValue, doF, numObs);
        pW_s = 10.^(-logp_s);
        pW_s(isinf(logp_s)) = realmin;

        % Expand to volume
        vol_W_s     = fullvol(W_s,     mask);
        vol_pW_s    = fullvol(pW_s,    mask);
        vol_logp_s  = fullvol(logp_s,  mask);

        % Bonferroni
        vol_pW_s_Bon          = min(max(vol_pW_s * correction_factor, min_pval_threshold), 1);
        vol_pW_s_Bon(mask<.5) = 1;
        vol_logp_s_Bon        = -log10(vol_pW_s_Bon);

        % FDR correction
        pvals_masked = pW_s;
        [h_FDR, crit_p] = fdr_bh(pvals_masked, 0.05, 'pdep');
        vol_pW_s_FDR = ones(size(mask));
        vol_pW_s_FDR(mask >= 0.5) = double(~h_FDR);  % 1 = not significant, 0 = significant
        % For visualization: use uncorrected -log10(p) where FDR is significant, 0 elsewhere
        vol_logp_s_FDR = vol_logp_s;  % Start with uncorrected log p-values
        vol_logp_s_FDR(vol_pW_s_FDR > 0.5) = 0;  % Set non-significant regions to 0
    end

    % ===========================
    % COMBINED LINEAR + SPLINE
    % ===========================
    if has_linear && has_splines
        fprintf('Performing combined effect analysis for: %s (Modality: %s)\n', xname, modality);

        statmat_combined = beta_hat([idx_splines, idx_linear, idx_intercept], :);
        % preallocate
        volmat_combined = zeros(size(bfmat,1), size(beta_hat,2), 'like', beta_hat);
        semat_combined  = zeros(size(bfmat,1), size(beta_hat,2), 'like', beta_hat);

        for xi = 1:size(bfmat, 1)
            wvec       = bfmat(xi, :)';
            % spline + linear + intercept (weights: [w; 1; 1])
            valvec_bf  = sum(statmat_combined .* [wvec; 1; 1], 1);
            sevec_bf   = sqrt( ...
                               pagemtimes( ...
                                  pagemtimes([wvec; 1; 1]', coeffCovar([idx_splines idx_linear idx_intercept], ...
                                                                       [idx_splines idx_linear idx_intercept], :)), ...
                                  [wvec; 1; 1] ...
                               ) ...
                             );
            volmat_combined(xi, :) = valvec_bf;
            semat_combined(xi, :)  = sevec_bf;
        end

        % Joint Wald (linear + all spline coeffs)
        comb_idx = [idx_linear, idx_splines];
        L = zeros(size(beta_hat, 1));
        L(comb_idx, comb_idx) = eye(length(comb_idx));

        hypValue = 0; doF = false; numObs = size(beta_hat, 2);
        [W_ls, logp_ls, ~, ~, dfW_ls, LW_ls] = FEMA_WaldTest(L, beta_hat, coeffCovar, hypValue, doF, numObs);
        pW_ls = 10.^(-logp_ls);
        pW_ls(isinf(logp_ls)) = realmin;

        % Expand to volume
        vol_W_ls     = fullvol(W_ls,     mask);
        vol_pW_ls    = fullvol(pW_ls,    mask);
        vol_logp_ls  = fullvol(logp_ls,  mask);

        % Bonferroni
        vol_pW_ls_Bon          = min(max(vol_pW_ls * correction_factor, min_pval_threshold), 1);
        vol_pW_ls_Bon(mask<.5) = 1;
        vol_logp_ls_Bon        = -log10(vol_pW_ls_Bon);

        % FDR correction
        pvals_masked = pW_ls;
        [h_FDR, crit_p] = fdr_bh(pvals_masked, 0.05, 'pdep');
        vol_pW_ls_FDR = ones(size(mask));
        vol_pW_ls_FDR(mask >= 0.5) = double(~h_FDR);  % 1 = not significant, 0 = significant
        % For visualization: use uncorrected -log10(p) where FDR is significant, 0 elsewhere
        vol_logp_ls_FDR = vol_logp_ls;  % Start with uncorrected log p-values
        vol_logp_ls_FDR(vol_pW_ls_FDR > 0.5) = 0;  % Set non-significant regions to 0
    end

    % ===========================
    % SAVE RESULTS
    % ===========================
    fname_out = sprintf('%s/%s_Wald_splines_CI_output_%s.mat', dirname_out, modality, xname);
    suffix    = ['_' xname];

    save_vars = {};

    if has_linear
        eval(['volmat_linear' suffix '   = volmat_linear;']);
        eval(['semat_linear'  suffix '   = semat_linear;']);
        eval(['W_l'           suffix '   = W_l;']);
        eval(['vol_W_l'       suffix '   = vol_W_l;']);
        eval(['pW_l'          suffix '   = pW_l;']);
        eval(['vol_pW_l'      suffix '   = vol_pW_l;']);
        eval(['logp_l'        suffix '   = logp_l;']);
        eval(['vol_logp_l'    suffix '   = vol_logp_l;']);
        eval(['dfW_l'         suffix '   = dfW_l;']);
        eval(['LW_l'          suffix '   = LW_l;']);
        eval(['vol_pW_l_Bon'  suffix '   = vol_pW_l_Bon;']);
        eval(['vol_logp_l_Bon' suffix '  = vol_logp_l_Bon;']);
        eval(['vol_pW_l_FDR'  suffix '   = vol_pW_l_FDR;']);
        eval(['vol_logp_l_FDR' suffix '  = vol_logp_l_FDR;']);
        save_vars = [save_vars, 'volmat_linear', 'semat_linear', 'W_l', 'vol_W_l', ...
                                'pW_l', 'vol_pW_l', 'logp_l', 'vol_logp_l', ...
                                'dfW_l', 'LW_l', 'vol_pW_l_Bon', 'vol_logp_l_Bon', ...
                                'vol_pW_l_FDR', 'vol_logp_l_FDR'];
    end

    if has_splines
        eval(['volmat_spline' suffix '   = volmat_spline;']);
        eval(['semat_spline'  suffix '   = semat_spline;']);
        eval(['W_s'           suffix '   = W_s;']);
        eval(['vol_W_s'       suffix '   = vol_W_s;']);
        eval(['pW_s'          suffix '   = pW_s;']);
        eval(['vol_pW_s'      suffix '   = vol_pW_s;']);
        eval(['logp_s'        suffix '   = logp_s;']);
        eval(['vol_logp_s'    suffix '   = vol_logp_s;']);
        eval(['dfW_s'         suffix '   = dfW_s;']);
        eval(['LW_s'          suffix '   = LW_s;']);
        eval(['vol_pW_s_Bon'  suffix '   = vol_pW_s_Bon;']);
        eval(['vol_logp_s_Bon' suffix '  = vol_logp_s_Bon;']);
        eval(['vol_pW_s_FDR'  suffix '   = vol_pW_s_FDR;']);
        eval(['vol_logp_s_FDR' suffix '  = vol_logp_s_FDR;']);
        save_vars = [save_vars, 'volmat_spline', 'semat_spline', 'W_s', 'vol_W_s', ...
                                'pW_s', 'vol_pW_s', 'logp_s', 'vol_logp_s', ...
                                'dfW_s', 'LW_s', 'vol_pW_s_Bon', 'vol_logp_s_Bon', ...
                                'vol_pW_s_FDR', 'vol_logp_s_FDR'];
    end

    if has_linear && has_splines
        eval(['volmat_combined' suffix '   = volmat_combined;']);
        eval(['semat_combined'  suffix '   = semat_combined;']);
        eval(['W_ls'            suffix '   = W_ls;']);
        eval(['vol_W_ls'        suffix '   = vol_W_ls;']);
        eval(['pW_ls'           suffix '   = pW_ls;']);
        eval(['vol_pW_ls'       suffix '   = vol_pW_ls;']);
        eval(['logp_ls'         suffix '   = logp_ls;']);
        eval(['vol_logp_ls'     suffix '   = vol_logp_ls;']);
        eval(['dfW_ls'          suffix '   = dfW_ls;']);
        eval(['LW_ls'           suffix '   = LW_ls;']);
        eval(['vol_pW_ls_Bon'   suffix '   = vol_pW_ls_Bon;']);
        eval(['vol_logp_ls_Bon' suffix '   = vol_logp_ls_Bon;']);
        eval(['vol_pW_ls_FDR'   suffix '   = vol_pW_ls_FDR;']);
        eval(['vol_logp_ls_FDR' suffix '   = vol_logp_ls_FDR;']);
        save_vars = [save_vars, 'volmat_combined', 'semat_combined', 'W_ls', 'vol_W_ls', ...
                                'pW_ls', 'vol_pW_ls', 'logp_ls', 'vol_logp_ls', ...
                                'dfW_ls', 'LW_ls', 'vol_pW_ls_Bon', 'vol_logp_ls_Bon', ...
                                'vol_pW_ls_FDR', 'vol_logp_ls_FDR'];
    end

    save(fname_out, save_vars{:}, '-v7.3');

    % Cleanup
    clear([save_vars{:} suffix]);
    clearvars -regexp '^volmat_' '^semat_' '^W_' '^pW_' '^logp_' '^LW_';
end