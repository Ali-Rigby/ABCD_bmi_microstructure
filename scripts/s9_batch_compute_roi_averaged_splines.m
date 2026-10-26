% ROI-averaged splines — unified BMI imaging pipeline (perc_median / bmip).

clear all; close all;

paths = get_project_paths();
addpath(paths.matlab_utils);
pheno_xname = paths.pheno_xname;
fprintf('s9 (pheno): BMI_PHENO_VARIABLE=%s\n', paths.pheno_variable);

% Define paths and file structure
release = "6.0";
dirname_outroot = paths.results_dir;
dirname_plot = fullfile(paths.results_dir, 'roi_data_for_plots');
fprintf('s9: results_dir = %s\n', dirname_outroot);
fprintf('s9: datafiles_dir = %s\n', paths.datafiles);
dirname_atlas = sprintf('/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0/concatenated/imaging/voxelwise/roi');

% Get design-matrix result folders only (exclude s9 plot output and other stray dirs)
designmat_dirs = dir(dirname_outroot);
designmat_dirs = designmat_dirs([designmat_dirs(:).isdir]);
names = {designmat_dirs(:).name};
designmat_dirs = designmat_dirs(startsWith(names, 'designmat'));

% Define model variables and effect types
variables = {paths.pheno_xname};
effect_types = {'splines'};

% Define atlas and ROI information
atlas_mod = {'aseg', 'fiber', 'thalamus_abcd3_cor10', 'pauli_subcortnuclei_abcd3_cor10'};
roi_mod = {
    {'Left-Thalamus-Proper', 'Right-Thalamus-Proper', 'Left-Caudate', 'Right-Caudate', 'Left-Putamen', 'Right-Putamen', ...
    'Left-Pallidum', 'Right-Pallidum', 'Left-Accumbens-area', 'Right-Accumbens-area', 'Left-Amygdala', 'Right-Amygdala', ...
    'Left-VentralDC', 'Right-VentralDC', 'Left-Cerebral-White-Matter', 'Right-Cerebral-White-Matter', 'Left-Cerebral-Cortex', ...
    'Right-Cerebral-Cortex', 'Left-Hippocampus', 'Right-Hippocampus', 'CC_Anterior'},  % for aseg 
    {'R_Fx', 'L_Fx', 'R_CgC', 'L_CgC', 'R_CgH', 'L_CgH', 'R_CST', 'L_CST', ...
    'R_ATR', 'L_ATR', 'R_Unc', 'L_Unc', 'R_ILF', 'L_ILF', 'R_IFO', 'L_IFO', ...
    'Fmaj', 'Fmin', 'CC', 'R_SLF', 'L_SLF', 'R_tSLF', 'L_tSLF', 'R_pSLF', 'L_pSLF', ...
    'R_SCS', 'L_SCS', 'R_fSCS', 'L_fSCS', 'R_pSCS', 'L_pSCS', 'R_SIFC', 'L_SIFC', ...
    'R_IFSFC', 'L_IFSFC'}, % for fiber
    {'LH-Pulvinar', 'LH-Anterior', 'LH-Medio_Dorsal', 'LH-Ventral_Latero_Dorsal', ... 
    'LH-Central_Lateral-Lateral_Posterior-Medial_Pulvinar', 'LH-Ventral_Anterior', 'LH-Ventral_Latero_Ventral', ... 
    'RH-Pulvinar', 'RH-Anterior', 'RH-Medio_Dorsal', 'RH-Ventral_Latero_Dorsal', 'RH-Central_Lateral-Lateral_Posterior-Medial_Pulvinar', ...
    'RH-Ventral_Anterior', 'RH-Ventral_Latero_Ventral'}, % for thalamus_abcd3_cor10
    {'Putamen','Caudate','Nucleus Accumbens','Extended Amygdala','Globus Pallidus external','Globus Pallidus internal', ...
    'Substantia Nigra pars compacta','Red Nucleus','Substantia Nigra pars reticulata','Parabrachial Pigmented Nucleus', ...
    'Ventral Tegmental Area','Ventral Pallidum','Habenular Nuclei','Hypothalamus','Mammillary Nucleus','Subthalamic Nucleus'} % for pauli_subcortnuclei_abcd3_cor10 
};

modlist = {'rni'};
datatype = 'voxel'; % must match s6
abcd_sync_path = fema_default_abcd_sync_path();

% Loop through each design matrix directory
for d = 1:length(designmat_dirs)
    dirname_desmat = designmat_dirs(d).name;
    dirname_out = fullfile(dirname_outroot, dirname_desmat);
    dirname_plot_out = fullfile(dirname_plot, dirname_desmat);

    % Print which design matrix folder is being processed
    fprintf('\n========================================\n');
    fprintf('Processing design matrix folder: %s\n', dirname_desmat);
    fprintf('========================================\n');

    if ~exist(dirname_plot_out, 'dir')
        mkdir(dirname_plot_out);
    end

    % Loop through each modality
    for mod_idx = 1:length(modlist)
        modality = modlist{mod_idx};
        
        fprintf('\n  Processing modality: %s\n', modality);
        
        dirname_imaging = fema_imaging_dirname(datatype, modality, abcd_sync_path);
        est_file = fema_results_path(dirname_out, datatype, modality);
        if ~isfile(est_file)
            warning('FEMA results file not found for %s - %s (%s), skipping...', ...
                dirname_desmat, modality, est_file);
            continue;
        end
        R = fema_load_results('ResultsFile', est_file, ...
            'Datatype', datatype, 'FstemImaging', modality, 'DirnameImaging', dirname_imaging);
        beta_hat = R.beta_hat;
        coeffCovar = R.coeffCovar;
        mask = R.mask;
        colnames_model = R.colnames_model;
        
        % 101-row evaluation grid (NOT *_interp.txt, which is per-subject BFs for s4)
        tbl_bf_pheno = readtable(fullfile(paths.datafiles, paths.pheno_bf_table));

        % Initialize tables to store summary results
        summary_table_data_Bon = cell(0, 5);
        summary_table_data_FDR = cell(0, 5);
 
% Loop over each variable (e.g., bmi, age, pds)
for var_idx = 1:length(variables)
    xname = variables{var_idx};
    
    % Define xvec based on variable
    switch xname
        case pheno_xname
            % Match s3 pheno grid (bmip: 0:100; perc_median: pheno_eval_grid)
            xvec = pheno_bf_xvec(tbl_bf_pheno, paths);
            fprintf('  %s xvec: n=%d, [%g : %g : %g]\n', xname, numel(xvec), xvec(1), xvec(2)-xvec(1), xvec(end));
        case 'age'
            xvec = (8.3:0.1:17.8);
        case 'pds'
            xvec = linspace(1, 5, 101);
    end

    % Load the corresponding results file for each variable
    fname_results = sprintf('%s/%s_Wald_splines_CI_output_%s.mat', dirname_out, modality, xname);
    if exist(fname_results, 'file')
        load(fname_results);
    else
        warning('Results file for %s not found, skipping...', xname);
        continue;
    end

    % Loop over each effect type: linear, splines
    for eff_idx = 1:length(effect_types)
        effect_type = effect_types{eff_idx};

        % Identify and load matrices based on effect type
        switch effect_type
            case 'linear'
                if exist('volmat_linear', 'var') && exist('semat_linear', 'var') && exist('vol_pW_l', 'var') && exist('vol_pW_l_FDR', 'var')
                    effect_volmat = volmat_linear; 
                    effect_semat = semat_linear; 
                    effect_vol_pW = vol_pW_l;
                    effect_vol_pW_FDR = vol_pW_l_FDR;
                else
                    continue;
                end
            case 'splines'
                if exist('volmat_spline', 'var') && exist('semat_spline', 'var') && exist('vol_pW_s', 'var') && exist('vol_pW_s_FDR', 'var')
                    effect_volmat = volmat_spline; 
                    effect_semat = semat_spline; 
                    effect_vol_pW = vol_pW_s;
                    effect_vol_pW_FDR = vol_pW_s_FDR;
                else
                    continue;
                end
        end

        % Define significance threshold and binary masks
        % Bonferroni correction
        pthresh = 0.05 / numel(find(mask >= 0.5));
        effect_vol_pW_bin_Bon = effect_vol_pW < pthresh;
        % FDR correction (already binary: 0=significant, 1=not significant)
        effect_vol_pW_bin_FDR = (effect_vol_pW_FDR == 0);
        
        % === COMPUTE DERIVATIVES FOR SPLINES ===
        if strcmp(effect_type, 'splines')
            fprintf('  Computing spline derivatives for %s...\n', xname);
            
            % Find spline coefficient indices based on variable name
            switch xname
                case pheno_xname
                    pat = sprintf('^%s_\\d+$', paths.pheno_bf_grep);
                    Svar_idx = find(find_cell(regexp(colnames_model, pat, 'match')));
                case 'age'
                    Svar_idx = find(find_cell(regexp(colnames_model,'^bf_age_\d+$','match')));
                case 'pds'
                    Svar_idx = find(find_cell(regexp(colnames_model,'^bf_pds_\d+$','match')));
            end
            intercept_idx = find(strcmp(colnames_model, 'intercept'));
            
            if ~isempty(Svar_idx)
                % Compute derivatives via FEMA_convert_splinevols
                try
                    switch xname
                        case pheno_xname
                            [~, ~, ~, volmat_dbf] = FEMA_convert_splinevols(tbl_bf_pheno, xvec, beta_hat, Svar_idx, intercept_idx, mask, coeffCovar);
                        case 'age'
                            [~, ~, ~, volmat_dbf] = FEMA_convert_splinevols(tbl_bf_age, xvec, beta_hat, Svar_idx, intercept_idx, mask, coeffCovar);
                        case 'pds'
                            [~, ~, ~, volmat_dbf] = FEMA_convert_splinevols(tbl_bf_pds, xvec, beta_hat, Svar_idx, intercept_idx, mask, coeffCovar);
                    end
                    
                    % Derivative mean matrix (xvec points x voxels)
                    effect_volmat_derivs = volmat_dbf;
                    
                    % Compute derivative standard errors manually
                    switch xname
                        case pheno_xname
                            bfmat = table2array(tbl_bf_pheno);
                        case 'age'
                            bfmat = table2array(tbl_bf_age);
                        case 'pds'
                            bfmat = table2array(tbl_bf_pds);
                    end
                    
                    nBasis = size(bfmat, 2);
                    if nBasis ~= numel(Svar_idx)
                        warning('Mismatch: table has %d basis columns, Svar_idx has %d betas. Derivative SEs may be invalid.', nBasis, numel(Svar_idx));
                        minBasis = min(nBasis, numel(Svar_idx));
                        bfmat = bfmat(:, 1:minBasis);
                    end
                    
                    % Gradient calculation (same as FEMA_convert_splinevols)
                    [~, dbfmat] = gradient(bfmat);
                    dbfmat = dbfmat / (xvec(2) - xvec(1));
                    nXvec = numel(xvec);
                    nVox = size(beta_hat, 2);
                    effect_semat_derivs = zeros(nXvec, nVox, 'single');
                    
                    % Vectorized covariance slice for selected betas
                    covBasis = coeffCovar(Svar_idx(1:size(bfmat,2)), Svar_idx(1:size(bfmat,2)), :);
                    for xi = 1:nXvec
                        wvec = dbfmat(xi, :)';
                        var_vec = pagemtimes(pagemtimes(wvec', covBasis), wvec);
                        se_vec = sqrt(squeeze(var_vec));
                        se_vec(~isnan(se_vec) & se_vec < 0 & se_vec > -1e-12) = 0;
                        effect_semat_derivs(xi, :) = se_vec;
                    end
                    
                    derivatives_computed = true;
                catch ME
                    warning('FEMA_convert_splinevols failed: %s', ME.message);
                    derivatives_computed = false;
                end
            else
                warning('No spline coefficients found for %s, skipping derivatives', xname);
                derivatives_computed = false;
            end
        else
            derivatives_computed = false;
        end

        % Loop through atlases and their respective ROIs
        for modi = 1:length(atlas_mod)
            atlas_type = atlas_mod{modi};
            fname_atlas = sprintf('%s/%s.mat', dirname_atlas, atlas_type);
            atlas = load(fname_atlas);
            
            % Extract ROI names and probability volumes based on atlas type
            if strcmp(atlas_type, 'thalamus_abcd3_cor10')
                roinames = atlas.thalamus_lut.LabelName;
                volmat_prob = atlas.vol_labels_reg;
            elseif strcmp(atlas_type, 'pauli_subcortnuclei_abcd3_cor10')
                roinames = atlas.pauli_lut.Name;
                volmat_prob = atlas.vol_labels_reg;
            else
                roinames = atlas.roinames;
                volmat_prob = atlas.volmat_prob;
            end
            
            strroi = roi_mod{modi};

            % Store bilateral ROI volumes for later processing
            bilateral_vols = struct();
            
            % Loop through each ROI
            for jj = 1:length(strroi)
                roi_name = strroi{jj};
                idx_roi = find(contains(roinames, roi_name));
                
                if isempty(idx_roi)
                    warning('ROI %s not found in atlas %s', roi_name, atlas_type);
                    continue;
                end

                % ROI-specific voxel processing with subsampling and masking
                vol_roi = subsample_volume(sum(volmat_prob(:,:,:,idx_roi) > 0.7, 4));
                vol_roi = imbinarize(vol_roi);
                
                % Store volume for bilateral processing
                if contains(roi_name, 'Left-') || contains(roi_name, 'LH-')
                    bilateral_name = strrep(strrep(roi_name, 'Left-', ''), 'LH-', '');
                    bilateral_name = strrep(bilateral_name, '-', '_');  % Replace hyphens with underscores for valid field names
                    bilateral_vols.(bilateral_name).left = vol_roi;
                    bilateral_vols.(bilateral_name).original_name = strrep(strrep(roi_name, 'Left-', ''), 'LH-', '');  % Store original name for output
                elseif contains(roi_name, 'Right-') || contains(roi_name, 'RH-')
                    bilateral_name = strrep(strrep(roi_name, 'Right-', ''), 'RH-', '');
                    bilateral_name = strrep(bilateral_name, '-', '_');  % Replace hyphens with underscores for valid field names
                    bilateral_vols.(bilateral_name).right = vol_roi;
                    bilateral_vols.(bilateral_name).original_name = strrep(strrep(roi_name, 'Right-', ''), 'RH-', '');  % Store original name for output
                elseif strcmp(atlas_type, 'fiber') && startsWith(roi_name, 'L_')
                    bilateral_name = matlab.lang.makeValidName(roi_name(3:end));
                    bilateral_vols.(bilateral_name).left = vol_roi;
                    bilateral_vols.(bilateral_name).original_name = roi_name(3:end);
                elseif strcmp(atlas_type, 'fiber') && startsWith(roi_name, 'R_')
                    bilateral_name = matlab.lang.makeValidName(roi_name(3:end));
                    bilateral_vols.(bilateral_name).right = vol_roi;
                    bilateral_vols.(bilateral_name).original_name = roi_name(3:end);
                end

                % Process all voxels in the ROI
                vol_roi_all = vol_roi .* mask;
                vec_roi_all = vol_roi_all(mask >= 0.5);
                vox_roi_all = find(vec_roi_all);
                num_voxels_total = length(vox_roi_all);

                % Compute mask for significant voxels within the ROI - Bonferroni
                vol_roi_sig_Bon = effect_vol_pW_bin_Bon .* vol_roi;
                vec_roi_inc_Bon = vol_roi_sig_Bon(mask >= 0.5);
                vox_roi_inc_Bon = find(vec_roi_inc_Bon);
                num_voxels_sig_Bon = length(vox_roi_inc_Bon);
                
                % Compute mask for significant voxels within the ROI - FDR
                vol_roi_sig_FDR = effect_vol_pW_bin_FDR .* vol_roi;
                vec_roi_inc_FDR = vol_roi_sig_FDR(mask >= 0.5);
                vox_roi_inc_FDR = find(vec_roi_inc_FDR);
                num_voxels_sig_FDR = length(vox_roi_inc_FDR);
                
                % Calculate percentage of significant voxels
                if num_voxels_total > 0
                    percent_sig_Bon = (num_voxels_sig_Bon / num_voxels_total) * 100;
                    percent_sig_FDR = (num_voxels_sig_FDR / num_voxels_total) * 100;
                else
                    percent_sig_Bon = 0;
                    percent_sig_FDR = 0;
                end
                
                % Add results to the summary tables
                summary_table_data_Bon(end+1, :) = {atlas_type, roi_name, num_voxels_total, num_voxels_sig_Bon, percent_sig_Bon};
                summary_table_data_FDR(end+1, :) = {atlas_type, roi_name, num_voxels_total, num_voxels_sig_FDR, percent_sig_FDR};

                % === BONFERRONI-CORRECTED SIGNIFICANT VOXELS ===
                effect_volmat_sig_Bon = effect_volmat(:, vox_roi_inc_Bon);
                effect_semat_sig_Bon = effect_semat(:, vox_roi_inc_Bon);

                if ~isempty(vox_roi_inc_Bon)
                    spline_ci_sig_Bon = NaN([2, size(effect_volmat, 1), length(vox_roi_inc_Bon)]);
                    for ii = 1:length(vox_roi_inc_Bon)
                        spline_ci_sig_Bon(1, :, ii) = effect_volmat_sig_Bon(:, ii) - 1.96 * effect_semat_sig_Bon(:, ii);
                        spline_ci_sig_Bon(2, :, ii) = effect_volmat_sig_Bon(:, ii) + 1.96 * effect_semat_sig_Bon(:, ii);
                    end
                    spline_ci_mean_sig_Bon = mean(spline_ci_sig_Bon, 3);
                    effect_volmat_mean_sig_Bon = mean(effect_volmat_sig_Bon, 2);
                else
                    spline_ci_mean_sig_Bon = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_sig_Bon = NaN(size(effect_volmat, 1), 1);
                end
                
                % === FDR-CORRECTED SIGNIFICANT VOXELS ===
                effect_volmat_sig_FDR = effect_volmat(:, vox_roi_inc_FDR);
                effect_semat_sig_FDR = effect_semat(:, vox_roi_inc_FDR);

                if ~isempty(vox_roi_inc_FDR)
                    spline_ci_sig_FDR = NaN([2, size(effect_volmat, 1), length(vox_roi_inc_FDR)]);
                    for ii = 1:length(vox_roi_inc_FDR)
                        spline_ci_sig_FDR(1, :, ii) = effect_volmat_sig_FDR(:, ii) - 1.96 * effect_semat_sig_FDR(:, ii);
                        spline_ci_sig_FDR(2, :, ii) = effect_volmat_sig_FDR(:, ii) + 1.96 * effect_semat_sig_FDR(:, ii);
                    end
                    spline_ci_mean_sig_FDR = mean(spline_ci_sig_FDR, 3);
                    effect_volmat_mean_sig_FDR = mean(effect_volmat_sig_FDR, 2);
                else
                    spline_ci_mean_sig_FDR = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_sig_FDR = NaN(size(effect_volmat, 1), 1);
                end

                effect_volmat_all = effect_volmat(:, vox_roi_all);
                effect_semat_all = effect_semat(:, vox_roi_all);

                if ~isempty(vox_roi_all)
                    spline_ci_all = NaN([2, size(effect_volmat, 1), length(vox_roi_all)]);
                    for ii = 1:length(vox_roi_all)
                        spline_ci_all(1, :, ii) = effect_volmat_all(:, ii) - 1.96 * effect_semat_all(:, ii);
                        spline_ci_all(2, :, ii) = effect_volmat_all(:, ii) + 1.96 * effect_semat_all(:, ii);
                    end
                    spline_ci_mean_all = mean(spline_ci_all, 3);
                    effect_volmat_mean_all = mean(effect_volmat_all, 2);
                else
                    spline_ci_mean_all = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_all = NaN(size(effect_volmat, 1), 1);
                end

                % Determine if this ROI is lateralized or non-lateralized
                is_lateralized = contains(roi_name, 'Left-') || contains(roi_name, 'Right-') || contains(roi_name, 'LH-') || contains(roi_name, 'RH-') || contains(roi_name, 'L_') || contains(roi_name, 'R_');
                if is_lateralized
                    laterality_dir = 'lateralized';
                else
                    laterality_dir = 'non_lateralized';
                end

                % Save Bonferroni-corrected significant results
                dirname_plot_Bon = fullfile(dirname_plot_out, 'betas', laterality_dir, 'bonferroni');
                if ~exist(dirname_plot_Bon, 'dir'), mkdir(dirname_plot_Bon); end
                txt_file_out_sig_Bon = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_Bon, modality, effect_type, xname, roi_name);
                fid_sig_Bon = fopen(txt_file_out_sig_Bon, 'w');
                if fid_sig_Bon == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_sig_Bon);
                end
                fprintf(fid_sig_Bon, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_sig_Bon, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_sig_Bon(xi), spline_ci_mean_sig_Bon(1, xi), spline_ci_mean_sig_Bon(2, xi));
                end
                fclose(fid_sig_Bon);
                
                % Save FDR-corrected significant results
                dirname_plot_FDR = fullfile(dirname_plot_out, 'betas', laterality_dir, 'fdr');
                if ~exist(dirname_plot_FDR, 'dir'), mkdir(dirname_plot_FDR); end
                txt_file_out_sig_FDR = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_FDR, modality, effect_type, xname, roi_name);
                fid_sig_FDR = fopen(txt_file_out_sig_FDR, 'w');
                if fid_sig_FDR == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_sig_FDR);
                end
                fprintf(fid_sig_FDR, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_sig_FDR, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_sig_FDR(xi), spline_ci_mean_sig_FDR(1, xi), spline_ci_mean_sig_FDR(2, xi));
                end
                fclose(fid_sig_FDR);

                % Save results for all voxels in the ROI
                dirname_plot_all = fullfile(dirname_plot_out, 'betas', laterality_dir, 'all_voxels');
                if ~exist(dirname_plot_all, 'dir'), mkdir(dirname_plot_all); end
                txt_file_out_all = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_all, modality, effect_type, xname, roi_name);
                fid_all = fopen(txt_file_out_all, 'w');
                if fid_all == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_all);
                end
                fprintf(fid_all, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_all, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_all(xi), spline_ci_mean_all(1, xi), spline_ci_mean_all(2, xi));
                end
                fclose(fid_all);
                
                % === SAVE DERIVATIVES IF COMPUTED (LATERALIZED) ===
                if derivatives_computed
                    % DERIVATIVES: Bonferroni-corrected significant voxels
                    effect_derivs_sig_Bon = effect_volmat_derivs(:, vox_roi_inc_Bon);
                    effect_derivs_se_sig_Bon = effect_semat_derivs(:, vox_roi_inc_Bon);
                    
                    if ~isempty(vox_roi_inc_Bon)
                        derivs_ci_sig_Bon = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_inc_Bon)]);
                        for ii = 1:length(vox_roi_inc_Bon)
                            derivs_ci_sig_Bon(1, :, ii) = effect_derivs_sig_Bon(:, ii) - 1.96 * effect_derivs_se_sig_Bon(:, ii);
                            derivs_ci_sig_Bon(2, :, ii) = effect_derivs_sig_Bon(:, ii) + 1.96 * effect_derivs_se_sig_Bon(:, ii);
                        end
                        derivs_ci_mean_sig_Bon = mean(derivs_ci_sig_Bon, 3);
                        effect_derivs_mean_sig_Bon = mean(effect_derivs_sig_Bon, 2);
                    else
                        derivs_ci_mean_sig_Bon = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_sig_Bon = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save Bonferroni derivatives
                    dirname_derivs_Bon = fullfile(dirname_plot_out, 'derivatives', laterality_dir, 'bonferroni');
                    if ~exist(dirname_derivs_Bon, 'dir'), mkdir(dirname_derivs_Bon); end
                    txt_file_derivs_Bon = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_Bon, modality, effect_type, xname, roi_name);
                    fid_derivs_Bon = fopen(txt_file_derivs_Bon, 'w');
                    fprintf(fid_derivs_Bon, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_Bon, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_sig_Bon(xi), derivs_ci_mean_sig_Bon(1, xi), derivs_ci_mean_sig_Bon(2, xi));
                    end
                    fclose(fid_derivs_Bon);
                    
                    % DERIVATIVES: FDR-corrected significant voxels
                    effect_derivs_sig_FDR = effect_volmat_derivs(:, vox_roi_inc_FDR);
                    effect_derivs_se_sig_FDR = effect_semat_derivs(:, vox_roi_inc_FDR);
                    
                    if ~isempty(vox_roi_inc_FDR)
                        derivs_ci_sig_FDR = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_inc_FDR)]);
                        for ii = 1:length(vox_roi_inc_FDR)
                            derivs_ci_sig_FDR(1, :, ii) = effect_derivs_sig_FDR(:, ii) - 1.96 * effect_derivs_se_sig_FDR(:, ii);
                            derivs_ci_sig_FDR(2, :, ii) = effect_derivs_sig_FDR(:, ii) + 1.96 * effect_derivs_se_sig_FDR(:, ii);
                        end
                        derivs_ci_mean_sig_FDR = mean(derivs_ci_sig_FDR, 3);
                        effect_derivs_mean_sig_FDR = mean(effect_derivs_sig_FDR, 2);
                    else
                        derivs_ci_mean_sig_FDR = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_sig_FDR = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save FDR derivatives
                    dirname_derivs_FDR = fullfile(dirname_plot_out, 'derivatives', laterality_dir, 'fdr');
                    if ~exist(dirname_derivs_FDR, 'dir'), mkdir(dirname_derivs_FDR); end
                    txt_file_derivs_FDR = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_FDR, modality, effect_type, xname, roi_name);
                    fid_derivs_FDR = fopen(txt_file_derivs_FDR, 'w');
                    fprintf(fid_derivs_FDR, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_FDR, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_sig_FDR(xi), derivs_ci_mean_sig_FDR(1, xi), derivs_ci_mean_sig_FDR(2, xi));
                    end
                    fclose(fid_derivs_FDR);
                    
                    % DERIVATIVES: All voxels in ROI
                    effect_derivs_all = effect_volmat_derivs(:, vox_roi_all);
                    effect_derivs_se_all = effect_semat_derivs(:, vox_roi_all);
                    
                    if ~isempty(vox_roi_all)
                        derivs_ci_all = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_all)]);
                        for ii = 1:length(vox_roi_all)
                            derivs_ci_all(1, :, ii) = effect_derivs_all(:, ii) - 1.96 * effect_derivs_se_all(:, ii);
                            derivs_ci_all(2, :, ii) = effect_derivs_all(:, ii) + 1.96 * effect_derivs_se_all(:, ii);
                        end
                        derivs_ci_mean_all = mean(derivs_ci_all, 3);
                        effect_derivs_mean_all = mean(effect_derivs_all, 2);
                    else
                        derivs_ci_mean_all = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_all = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save all voxels derivatives
                    dirname_derivs_all = fullfile(dirname_plot_out, 'derivatives', laterality_dir, 'all_voxels');
                    if ~exist(dirname_derivs_all, 'dir'), mkdir(dirname_derivs_all); end
                    txt_file_derivs_all = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_all, modality, effect_type, xname, roi_name);
                    fid_derivs_all = fopen(txt_file_derivs_all, 'w');
                    fprintf(fid_derivs_all, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_all, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_all(xi), derivs_ci_mean_all(1, xi), derivs_ci_mean_all(2, xi));
                    end
                    fclose(fid_derivs_all);
                    
                    clear effect_derivs_sig_Bon effect_derivs_se_sig_Bon effect_derivs_sig_FDR effect_derivs_se_sig_FDR effect_derivs_all effect_derivs_se_all;
                    clear derivs_ci_sig_Bon derivs_ci_sig_FDR derivs_ci_all;
                end
                
                clear vol_roi vol_roi_sig_Bon vol_roi_sig_FDR vol_roi_all vec_roi_inc_Bon vec_roi_inc_FDR vec_roi_all vox_roi_inc_Bon vox_roi_inc_FDR vox_roi_all ...
                      effect_volmat_sig_Bon effect_semat_sig_Bon effect_volmat_sig_FDR effect_semat_sig_FDR effect_volmat_all effect_semat_all ...
                      spline_ci_sig_Bon spline_ci_sig_FDR spline_ci_all;
            end
            
            % === PROCESS BILATERAL STRUCTURES ===
            fprintf('\n  Processing bilateral structures for atlas: %s\n', atlas_type);
            bilateral_names = fieldnames(bilateral_vols);
            for bb = 1:length(bilateral_names)
                bilateral_name = bilateral_names{bb};
                
                % Check if both left and right exist
                if ~isfield(bilateral_vols.(bilateral_name), 'left') || ~isfield(bilateral_vols.(bilateral_name), 'right')
                    continue;
                end
                
                % Combine left and right volumes
                vol_roi = bilateral_vols.(bilateral_name).left | bilateral_vols.(bilateral_name).right;
                % Use original name (with hyphens) for output files
                roi_name = ['Bilateral-' bilateral_vols.(bilateral_name).original_name];
                
                fprintf('    Processing bilateral ROI: %s\n', roi_name);
                
                % Process all voxels in the bilateral ROI
                vol_roi_all = vol_roi .* mask;
                vec_roi_all = vol_roi_all(mask >= 0.5);
                vox_roi_all = find(vec_roi_all);
                num_voxels_total = length(vox_roi_all);

                % Compute mask for significant voxels within the ROI - Bonferroni
                vol_roi_sig_Bon = effect_vol_pW_bin_Bon .* vol_roi;
                vec_roi_inc_Bon = vol_roi_sig_Bon(mask >= 0.5);
                vox_roi_inc_Bon = find(vec_roi_inc_Bon);
                num_voxels_sig_Bon = length(vox_roi_inc_Bon);
                
                % Compute mask for significant voxels within the ROI - FDR
                vol_roi_sig_FDR = effect_vol_pW_bin_FDR .* vol_roi;
                vec_roi_inc_FDR = vol_roi_sig_FDR(mask >= 0.5);
                vox_roi_inc_FDR = find(vec_roi_inc_FDR);
                num_voxels_sig_FDR = length(vox_roi_inc_FDR);
                
                % Calculate percentage of significant voxels
                if num_voxels_total > 0
                    percent_sig_Bon = (num_voxels_sig_Bon / num_voxels_total) * 100;
                    percent_sig_FDR = (num_voxels_sig_FDR / num_voxels_total) * 100;
                else
                    percent_sig_Bon = 0;
                    percent_sig_FDR = 0;
                end
                
                % Add results to the summary tables
                summary_table_data_Bon(end+1, :) = {atlas_type, roi_name, num_voxels_total, num_voxels_sig_Bon, percent_sig_Bon};
                summary_table_data_FDR(end+1, :) = {atlas_type, roi_name, num_voxels_total, num_voxels_sig_FDR, percent_sig_FDR};

                % === BONFERRONI-CORRECTED SIGNIFICANT VOXELS ===
                effect_volmat_sig_Bon = effect_volmat(:, vox_roi_inc_Bon);
                effect_semat_sig_Bon = effect_semat(:, vox_roi_inc_Bon);

                if ~isempty(vox_roi_inc_Bon)
                    spline_ci_sig_Bon = NaN([2, size(effect_volmat, 1), length(vox_roi_inc_Bon)]);
                    for ii = 1:length(vox_roi_inc_Bon)
                        spline_ci_sig_Bon(1, :, ii) = effect_volmat_sig_Bon(:, ii) - 1.96 * effect_semat_sig_Bon(:, ii);
                        spline_ci_sig_Bon(2, :, ii) = effect_volmat_sig_Bon(:, ii) + 1.96 * effect_semat_sig_Bon(:, ii);
                    end
                    spline_ci_mean_sig_Bon = mean(spline_ci_sig_Bon, 3);
                    effect_volmat_mean_sig_Bon = mean(effect_volmat_sig_Bon, 2);
                else
                    spline_ci_mean_sig_Bon = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_sig_Bon = NaN(size(effect_volmat, 1), 1);
                end
                
                % === FDR-CORRECTED SIGNIFICANT VOXELS ===
                effect_volmat_sig_FDR = effect_volmat(:, vox_roi_inc_FDR);
                effect_semat_sig_FDR = effect_semat(:, vox_roi_inc_FDR);

                if ~isempty(vox_roi_inc_FDR)
                    spline_ci_sig_FDR = NaN([2, size(effect_volmat, 1), length(vox_roi_inc_FDR)]);
                    for ii = 1:length(vox_roi_inc_FDR)
                        spline_ci_sig_FDR(1, :, ii) = effect_volmat_sig_FDR(:, ii) - 1.96 * effect_semat_sig_FDR(:, ii);
                        spline_ci_sig_FDR(2, :, ii) = effect_volmat_sig_FDR(:, ii) + 1.96 * effect_semat_sig_FDR(:, ii);
                    end
                    spline_ci_mean_sig_FDR = mean(spline_ci_sig_FDR, 3);
                    effect_volmat_mean_sig_FDR = mean(effect_volmat_sig_FDR, 2);
                else
                    spline_ci_mean_sig_FDR = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_sig_FDR = NaN(size(effect_volmat, 1), 1);
                end

                effect_volmat_all = effect_volmat(:, vox_roi_all);
                effect_semat_all = effect_semat(:, vox_roi_all);

                if ~isempty(vox_roi_all)
                    spline_ci_all = NaN([2, size(effect_volmat, 1), length(vox_roi_all)]);
                    for ii = 1:length(vox_roi_all)
                        spline_ci_all(1, :, ii) = effect_volmat_all(:, ii) - 1.96 * effect_semat_all(:, ii);
                        spline_ci_all(2, :, ii) = effect_volmat_all(:, ii) + 1.96 * effect_semat_all(:, ii);
                    end
                    spline_ci_mean_all = mean(spline_ci_all, 3);
                    effect_volmat_mean_all = mean(effect_volmat_all, 2);
                else
                    spline_ci_mean_all = NaN(2, size(effect_volmat, 1));
                    effect_volmat_mean_all = NaN(size(effect_volmat, 1), 1);
                end

                % Save Bonferroni-corrected significant results
                dirname_plot_Bon = fullfile(dirname_plot_out, 'betas', 'bilateral', 'bonferroni');
                if ~exist(dirname_plot_Bon, 'dir'), mkdir(dirname_plot_Bon); end
                txt_file_out_sig_Bon = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_Bon, modality, effect_type, xname, roi_name);
                fid_sig_Bon = fopen(txt_file_out_sig_Bon, 'w');
                if fid_sig_Bon == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_sig_Bon);
                end
                fprintf(fid_sig_Bon, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_sig_Bon, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_sig_Bon(xi), spline_ci_mean_sig_Bon(1, xi), spline_ci_mean_sig_Bon(2, xi));
                end
                fclose(fid_sig_Bon);
                
                % Save FDR-corrected significant results
                dirname_plot_FDR = fullfile(dirname_plot_out, 'betas', 'bilateral', 'fdr');
                if ~exist(dirname_plot_FDR, 'dir'), mkdir(dirname_plot_FDR); end
                txt_file_out_sig_FDR = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_FDR, modality, effect_type, xname, roi_name);
                fid_sig_FDR = fopen(txt_file_out_sig_FDR, 'w');
                if fid_sig_FDR == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_sig_FDR);
                end
                fprintf(fid_sig_FDR, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_sig_FDR, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_sig_FDR(xi), spline_ci_mean_sig_FDR(1, xi), spline_ci_mean_sig_FDR(2, xi));
                end
                fclose(fid_sig_FDR);

                % Save results for all voxels in the ROI
                dirname_plot_all = fullfile(dirname_plot_out, 'betas', 'bilateral', 'all_voxels');
                if ~exist(dirname_plot_all, 'dir'), mkdir(dirname_plot_all); end
                txt_file_out_all = sprintf('%s/%s_%s_%s_CI_%s.txt', dirname_plot_all, modality, effect_type, xname, roi_name);
                fid_all = fopen(txt_file_out_all, 'w');
                if fid_all == -1
                    error('Error opening file %s for writing. Check directory permissions.', txt_file_out_all);
                end
                fprintf(fid_all, 'xvec,mean,CI_lower,CI_upper\n');
                for xi = 1:length(xvec)
                    fprintf(fid_all, '%f,%f,%f,%f\n', xvec(xi), effect_volmat_mean_all(xi), spline_ci_mean_all(1, xi), spline_ci_mean_all(2, xi));
                end
                fclose(fid_all);
                
                % === SAVE DERIVATIVES IF COMPUTED (BILATERAL) ===
                if derivatives_computed
                    % DERIVATIVES: Bonferroni-corrected significant voxels
                    effect_derivs_sig_Bon = effect_volmat_derivs(:, vox_roi_inc_Bon);
                    effect_derivs_se_sig_Bon = effect_semat_derivs(:, vox_roi_inc_Bon);
                    
                    if ~isempty(vox_roi_inc_Bon)
                        derivs_ci_sig_Bon = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_inc_Bon)]);
                        for ii = 1:length(vox_roi_inc_Bon)
                            derivs_ci_sig_Bon(1, :, ii) = effect_derivs_sig_Bon(:, ii) - 1.96 * effect_derivs_se_sig_Bon(:, ii);
                            derivs_ci_sig_Bon(2, :, ii) = effect_derivs_sig_Bon(:, ii) + 1.96 * effect_derivs_se_sig_Bon(:, ii);
                        end
                        derivs_ci_mean_sig_Bon = mean(derivs_ci_sig_Bon, 3);
                        effect_derivs_mean_sig_Bon = mean(effect_derivs_sig_Bon, 2);
                    else
                        derivs_ci_mean_sig_Bon = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_sig_Bon = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save Bonferroni derivatives
                    dirname_derivs_Bon = fullfile(dirname_plot_out, 'derivatives', 'bilateral', 'bonferroni');
                    if ~exist(dirname_derivs_Bon, 'dir'), mkdir(dirname_derivs_Bon); end
                    txt_file_derivs_Bon = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_Bon, modality, effect_type, xname, roi_name);
                    fid_derivs_Bon = fopen(txt_file_derivs_Bon, 'w');
                    fprintf(fid_derivs_Bon, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_Bon, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_sig_Bon(xi), derivs_ci_mean_sig_Bon(1, xi), derivs_ci_mean_sig_Bon(2, xi));
                    end
                    fclose(fid_derivs_Bon);
                    
                    % DERIVATIVES: FDR-corrected significant voxels
                    effect_derivs_sig_FDR = effect_volmat_derivs(:, vox_roi_inc_FDR);
                    effect_derivs_se_sig_FDR = effect_semat_derivs(:, vox_roi_inc_FDR);
                    
                    if ~isempty(vox_roi_inc_FDR)
                        derivs_ci_sig_FDR = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_inc_FDR)]);
                        for ii = 1:length(vox_roi_inc_FDR)
                            derivs_ci_sig_FDR(1, :, ii) = effect_derivs_sig_FDR(:, ii) - 1.96 * effect_derivs_se_sig_FDR(:, ii);
                            derivs_ci_sig_FDR(2, :, ii) = effect_derivs_sig_FDR(:, ii) + 1.96 * effect_derivs_se_sig_FDR(:, ii);
                        end
                        derivs_ci_mean_sig_FDR = mean(derivs_ci_sig_FDR, 3);
                        effect_derivs_mean_sig_FDR = mean(effect_derivs_sig_FDR, 2);
                    else
                        derivs_ci_mean_sig_FDR = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_sig_FDR = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save FDR derivatives
                    dirname_derivs_FDR = fullfile(dirname_plot_out, 'derivatives', 'bilateral', 'fdr');
                    if ~exist(dirname_derivs_FDR, 'dir'), mkdir(dirname_derivs_FDR); end
                    txt_file_derivs_FDR = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_FDR, modality, effect_type, xname, roi_name);
                    fid_derivs_FDR = fopen(txt_file_derivs_FDR, 'w');
                    fprintf(fid_derivs_FDR, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_FDR, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_sig_FDR(xi), derivs_ci_mean_sig_FDR(1, xi), derivs_ci_mean_sig_FDR(2, xi));
                    end
                    fclose(fid_derivs_FDR);
                    
                    % DERIVATIVES: All voxels in ROI
                    effect_derivs_all = effect_volmat_derivs(:, vox_roi_all);
                    effect_derivs_se_all = effect_semat_derivs(:, vox_roi_all);
                    
                    if ~isempty(vox_roi_all)
                        derivs_ci_all = NaN([2, size(effect_volmat_derivs, 1), length(vox_roi_all)]);
                        for ii = 1:length(vox_roi_all)
                            derivs_ci_all(1, :, ii) = effect_derivs_all(:, ii) - 1.96 * effect_derivs_se_all(:, ii);
                            derivs_ci_all(2, :, ii) = effect_derivs_all(:, ii) + 1.96 * effect_derivs_se_all(:, ii);
                        end
                        derivs_ci_mean_all = mean(derivs_ci_all, 3);
                        effect_derivs_mean_all = mean(effect_derivs_all, 2);
                    else
                        derivs_ci_mean_all = NaN(2, size(effect_volmat_derivs, 1));
                        effect_derivs_mean_all = NaN(size(effect_volmat_derivs, 1), 1);
                    end
                    
                    % Save all voxels derivatives
                    dirname_derivs_all = fullfile(dirname_plot_out, 'derivatives', 'bilateral', 'all_voxels');
                    if ~exist(dirname_derivs_all, 'dir'), mkdir(dirname_derivs_all); end
                    txt_file_derivs_all = sprintf('%s/%s_%s_%s_derivatives_CI_%s.txt', dirname_derivs_all, modality, effect_type, xname, roi_name);
                    fid_derivs_all = fopen(txt_file_derivs_all, 'w');
                    fprintf(fid_derivs_all, 'xvec,mean_derivative,CI_lower,CI_upper\n');
                    for xi = 1:length(xvec)
                        fprintf(fid_derivs_all, '%f,%f,%f,%f\n', xvec(xi), effect_derivs_mean_all(xi), derivs_ci_mean_all(1, xi), derivs_ci_mean_all(2, xi));
                    end
                    fclose(fid_derivs_all);
                    
                    clear effect_derivs_sig_Bon effect_derivs_se_sig_Bon effect_derivs_sig_FDR effect_derivs_se_sig_FDR effect_derivs_all effect_derivs_se_all;
                    clear derivs_ci_sig_Bon derivs_ci_sig_FDR derivs_ci_all;
                end
                
                clear vol_roi vol_roi_sig_Bon vol_roi_sig_FDR vol_roi_all vec_roi_inc_Bon vec_roi_inc_FDR vec_roi_all vox_roi_inc_Bon vox_roi_inc_FDR vox_roi_all ...
                      effect_volmat_sig_Bon effect_semat_sig_Bon effect_volmat_sig_FDR effect_semat_sig_FDR effect_volmat_all effect_semat_all ...
                      spline_ci_sig_Bon spline_ci_sig_FDR spline_ci_all;
            end
            
            clear fname_atlas atlas strroi bilateral_vols;
        end
        clear effect_volmat effect_semat effect_vol_pW effect_vol_pW_FDR effect_vol_pW_bin_Bon effect_vol_pW_bin_FDR;
    end
    clear xname xvec fname_results;
end

    % Display and save summary tables for this modality
    fprintf('\n--- ROI Significance Summary (Bonferroni) for %s - %s ---\n', dirname_desmat, modality);
    summary_table_Bon = cell2table(summary_table_data_Bon, 'VariableNames', {'Atlas', 'ROI', 'TotalVoxels', 'SignificantVoxels', 'PercentSignificant'});
    disp(summary_table_Bon);
    
    fprintf('\n--- ROI Significance Summary (FDR) for %s - %s ---\n', dirname_desmat, modality);
    summary_table_FDR = cell2table(summary_table_data_FDR, 'VariableNames', {'Atlas', 'ROI', 'TotalVoxels', 'SignificantVoxels', 'PercentSignificant'});
    disp(summary_table_FDR);
    
    % Save summary tables as CSV files
    csv_file_Bon = sprintf('%s/ROI_summary_Bonferroni_sig_voxels_%s.csv', dirname_plot_out, modality);
    writetable(summary_table_Bon, csv_file_Bon);
    fprintf('Saved Bonferroni summary to: %s\n', csv_file_Bon);
    
    csv_file_FDR = sprintf('%s/ROI_summary_FDR_sig_voxels_%s.csv', dirname_plot_out, modality);
    writetable(summary_table_FDR, csv_file_FDR);
    fprintf('Saved FDR summary to: %s\n', csv_file_FDR);

    % Clear summary tables for next modality
    clear summary_table_data_Bon summary_table_data_FDR summary_table_Bon summary_table_FDR;
    
    end  % End of modality loop

end  % End of design matrix loop