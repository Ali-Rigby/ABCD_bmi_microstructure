%% Build showVol volume bundles (Wald + spline derivatives) for screenshots.
% Unified BMI imaging pipeline (bmip / perc_median).
%
% Before running, source the matching env:
%   source env/run_bmip.csh        % or run_perc_median.csh
%
% Optional: set pheno_list below to {'bmip','perc_median'} to process both in one session.
% GUI-friendly: set pheno_list, results_run, and pheno_grid_step below (no .csh required).
% Optional env: BMI_IMAGING_RESULTS_RUN, BMI_PHENO_VARIABLE
%
% Outputs (under results dir): results_for_screenshots_voxel_<xname>/

clearvars -except pheno_list desmat_filter targetMod datatype summary_modality wald_corrections results_run pheno_grid_step;
close all;

%% User knobs
desmat_filter = '1m';   % '1f' | '1m' | '1' (1f+1m) | 'all'
targetMod     = 'rni';  % FEMA fstem_imaging / modality tag
datatype      = 'voxel';
pheno_list    = {'perc_median'};     % {} = use BMI_PHENO_VARIABLE from env; or {'bmip','perc_median'}
results_run   = 'results_2026_06_29'; % '' = env / auto date (see get_project_paths)
pheno_grid_step = 1;                % perc_median only; must match s3 (run_perc_median.csh uses 1)
wald_corrections = {'fdr'};  % subset of: 'uncorr', 'bonferroni', 'fdr'

%% Paths (this file lives in figures/matlab/)
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'scripts', 'matlab_add_paths.m'));

if isempty(pheno_list)
    pheno_env = strtrim(getenv('BMI_PHENO_VARIABLE'));
    if isempty(pheno_env)
        error(['Set BMI_PHENO_VARIABLE (source env/run_bmip.csh or run_perc_median.csh) ', ...
            'or set pheno_list = {''bmip''} / {''perc_median''} in this script.']);
    end
    pheno_list = {pheno_env};
end

%% Visual defaults (tune per modality / pheno as needed)
showPval = false;
cmap_linear = redblackblue;
cmap_splines = redblackblue;
cmap_wald = cmap_linear(200:end, :);
bg_vol = [];
CSFprob = [];
interpMethod = [];
atlasVersion = 'ABCD3';

pwaldlims_by_xname = struct('bmi', [0, 30], 'perc_median', [0, 30]);
limits_deriv_by_xname = struct('bmi', [-0.0009, 0.0009], 'perc_median', [-0.0009, 0.0009]);

for ip = 1:numel(pheno_list)
    pheno = strtrim(pheno_list{ip});
    setenv('BMI_PHENO_VARIABLE', pheno);
    setenv('BMI_IMAGING_ANALYSIS_DIR', pheno);
    if ~isempty(results_run)
        setenv('BMI_IMAGING_RESULTS_RUN', results_run);
    end
    if strcmp(pheno, 'perc_median') && ~isempty(pheno_grid_step)
        setenv('BMI_IMAGING_PHENO_GRID_STEP', num2str(pheno_grid_step));
    end

    paths = get_project_paths();
    addpath(genpath(paths.cmig_tools));
    xname = paths.pheno_xname;
    fprintf('\n=== showVol prep: %s (xname=%s) ===\n', pheno, xname);
    fprintf('Results dir: %s\n', paths.results_dir);

    if ~isfield(pwaldlims_by_xname, xname)
        error('Add pwaldlims for xname=%s in pwaldlims_by_xname.', xname);
    end
    if ~isfield(limits_deriv_by_xname, xname)
        error('Add derivative limits for xname=%s in limits_deriv_by_xname.', xname);
    end

    pwaldlims = pwaldlims_by_xname.(xname);
    pwaldlims_Bon = pwaldlims;
    pwaldlims_FDR = pwaldlims;
    limits_deriv = limits_deriv_by_xname.(xname);

    dirname_out = paths.results_dir;
    dirname_save = fullfile(dirname_out, sprintf('results_for_screenshots_%s_%s', datatype, xname));
    if ~exist(dirname_save, 'dir')
        mkdir(dirname_save);
    end

    designmat_names = list_designmat_dirs(dirname_out);
    if isempty(designmat_names)
        warning('No designmat* folders under %s', dirname_out);
        continue;
    end

    des_idx = select_designmat_indices(designmat_names, desmat_filter);
    if isempty(des_idx)
        error('No designmat folders match filter ''%s'' under %s', desmat_filter, dirname_out);
    end

    tbl_bf_pheno = readtable(fullfile(paths.datafiles, paths.pheno_bf_table), 'Delimiter', '\t');
    phenovec = pheno_bf_xvec(tbl_bf_pheno, paths);
    pheno_samples = 2:(numel(phenovec) - 1);

    abcd_sync_path = fema_default_abcd_sync_path();
    dirname_imaging = fema_imaging_dirname(datatype, targetMod, abcd_sync_path);

    for d = des_idx
        dname = designmat_names{d};
        dirname_results_d = fullfile(dirname_out, dname);
        desmat_tag = designmat_tag_from_name(dname);

        outdir_d = fullfile(dirname_save, dname);
        if ~exist(outdir_d, 'dir')
            mkdir(outdir_d);
        end

        est_file = fema_results_path(dirname_results_d, datatype, targetMod);
        R = fema_load_results('ResultsFile', est_file, ...
            'Datatype', datatype, 'FstemImaging', targetMod, ...
            'DirnameImaging', dirname_imaging);

        beta_hat = R.beta_hat;
        coeffCovar = R.coeffCovar;
        mask = R.mask;
        colnames_model = R.colnames_model;

        Spheno_idx = find(find_cell(regexp(colnames_model, ['^', paths.pheno_bf_grep], 'match')));
        intercept_idx = find(strcmp(colnames_model, 'intercept'));
        if isempty(Spheno_idx)
            error('No %s basis columns in %s', paths.pheno_bf_grep, est_file);
        end

        fprintf('Processing: %s | %s | %s\n', dname, datatype, targetMod);

        %% Wald maps
        vols_wald = process_wald_stats([], dirname_results_d, targetMod, xname, ...
            dname, pwaldlims, pwaldlims_Bon, pwaldlims_FDR, ...
            showPval, cmap_wald, bg_vol, CSFprob, interpMethod, atlasVersion, ...
            'Corrections', wald_corrections);
        labels_wald = build_wald_labels(desmat_tag, xname, wald_corrections, numel(vols_wald));

        if ~isempty(vols_wald)
            vols = vols_wald;
            labels = labels_wald;
            save(fullfile(outdir_d, sprintf('vols_summary_%s_wald.mat', targetMod)), ...
                'vols', 'labels', '-v7.3');
        end

        %% Spline derivative maps (mean derivative only; no CI bands)
        [vols_derivatives, ~, ~, ~] = FEMA_convert_splinevols( ...
            tbl_bf_pheno, phenovec, beta_hat, Spheno_idx, intercept_idx, mask, coeffCovar);

        vols_deriv = {};
        labels_deriv = {};
        for k = 1:numel(pheno_samples)
            slice = pheno_samples(k);
            pheno_val = phenovec(slice);
            statname_deriv = sprintf('beta_hat_derivative (%s = %s), %s', ...
                xname, format_pheno_value(pheno_val), dname);
            vols_deriv_converted = convert_and_append_vols([], vols_derivatives(:, :, :, slice), ...
                targetMod, statname_deriv, colnames_model(Spheno_idx), 1, limits_deriv, ...
                showPval, cmap_splines, bg_vol, CSFprob, interpMethod, atlasVersion);
            vols_deriv{end + 1} = vols_deriv_converted{1}; %#ok<AGROW>
            labels_deriv{end + 1} = sprintf('spline_deriv_%s_%s_%s', ...
                desmat_tag, xname, format_pheno_label(pheno_val)); %#ok<AGROW>
        end

        if ~isempty(vols_deriv)
            vols = vols_deriv;
            labels = labels_deriv;
            save(fullfile(outdir_d, sprintf('vols_summary_%s_derivatives.mat', targetMod)), ...
                'vols', 'labels', '-v7.3');
        end
    end

    fprintf('Saved screenshot bundles under: %s\n', dirname_save);
end

%% Local helpers
function names = list_designmat_dirs(dirname_out)
    d = dir(fullfile(dirname_out, 'designmat*'));
    d = d([d.isdir]);
    names = {d.name}';
    keep = false(size(names));
    for i = 1:numel(names)
        keep(i) = ~isempty(regexp(names{i}, '^designmat\d+[fm]_', 'once'));
    end
    names = names(keep);
    names = sort(names);
end

function idx = select_designmat_indices(designmat_names, desmat_filter)
    if strcmpi(strtrim(desmat_filter), 'all')
        idx = (1:numel(designmat_names))';
        return;
    end
    pattern = sprintf('designmat%s_', desmat_filter);
    idx = find(contains(designmat_names, pattern))';
end

function tag = designmat_tag_from_name(dname)
    tok = regexp(dname, '^designmat(\d+[fm])_', 'tokens', 'once');
    if isempty(tok)
        tag = 'unknown';
    else
        tag = tok{1};
    end
end

function s = format_pheno_value(v)
    if abs(v - round(v)) < 1e-9
        s = sprintf('%d', round(v));
    else
        s = sprintf('%.1f', v);
    end
end

function s = format_pheno_label(v)
    if v < 0
        s = sprintf('m%s', format_pheno_value(abs(v)));
    else
        s = format_pheno_label_pos(v);
    end
end

function s = format_pheno_label_pos(v)
    if abs(v - round(v)) < 1e-9
        s = sprintf('%d', round(v));
    else
        s = strrep(sprintf('%.1f', v), '.', 'p');
    end
end

function labels = build_wald_labels(desmat_tag, xname, corrections, nvols)
    stat_suffix = {'l', 's', 'ls'};
    labels = {};
    for i = 1:numel(corrections)
        corr = lower(strtrim(corrections{i}));
        for s = 1:numel(stat_suffix)
            labels{end + 1} = sprintf('wald_%s_%s_%s_%s', desmat_tag, xname, corr, stat_suffix{s}); %#ok<AGROW>
        end
    end
    if numel(labels) > nvols
        labels = labels(1:nvols);
    elseif numel(labels) < nvols
        for k = (numel(labels) + 1):nvols
            labels{k} = sprintf('wald_%s_%s_%02d', desmat_tag, xname, k);
        end
    end
end
