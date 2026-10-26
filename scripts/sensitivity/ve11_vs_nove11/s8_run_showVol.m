% Extended showVol for BMI imaging: Wald chi2 + -log10(p), pheno/age splines & derivatives.
%
% Keeps s8_run_showVol_pheno.m unchanged. Voxel maps only (no vertex branch).
%
% Path selection: same as s8 — set analysis_profile below or source env first.
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(this_dir)));  % .../scripts/sensitivity/ve11_vs_nove11 -> repo
run(fullfile(repo_root, 'scripts', 'matlab_add_paths.m'));

%% Analysis profile — edit these two lines in the GUI (non-empty overrides stale env)
analysis_profile = 've11_nove11_noVE11_perc_median';   % or ve11_nove11_noVE11_bmip | bmip_no_VE11 | r7_htTol5cm_VE11_bmip
results_run      = 'results_2026_10_02'; % ve11_vs_noVE11 FEMA/s7 run on ceph

if isempty(analysis_profile)
    analysis_profile = getenv('BMI_S8_ANALYSIS_PROFILE');
    if isempty(analysis_profile)
        analysis_profile = 've11_nove11_VE11_bmip';
    end
end
if isempty(results_run)
    results_run = getenv('BMI_IMAGING_RESULTS_RUN');
    if isempty(results_run)
        results_run = 'results_2026_09_02';
    end
end

%% Extended visualization knobs
% Predictors to visualize (default: active pheno + age). Must match s7 Wald stems.
% Leave empty ({}) to auto-fill from the analysis profile (bmip -> 'bmi', perc_median -> 'perc_median').
% If you hardcode 'bmi' here while using a perc_median profile, pheno maps will be skipped.
effect_vars = {};  % e.g. {'perc_median', 'age'} or {'bmi', 'age'}; {} = auto
wald_corrections = {'uncorr'};  % subset of: 'uncorr', 'bonferroni', 'fdr'
include_wald_stat = true;       % vol_W_* (chi-square)
include_logp = true;            % vol_logp_* (uses wald_corrections when true)
include_spline_slices = true;   % beta_hat slices along BF grid
include_derivatives = false;     % d(beta)/dx spline derivatives
include_sex_tstat =false;       % sexF t/z from FEMA zmat -> fullvol (combined-sex models)
vol_display_mode = 'color';    % 'both' | 'color' | 'gray' — applies to every map appended

% Spline / derivative slice selection (separate specs per predictor + map type).
% Modes: 'default' | 'indices' | 'values' | 'step' | 'decade' | 'all_interior'
% Age slices are profile-specific (see s8_init_slice_specs):
%   bmip / perc_median / scanner / ve11_nove11_* -> 8.5, 9:17, 17.8 on 8.3:0.1:17.8
%   bmip_no_VE11 -> 9:17 + 17.8 on 9.0:0.1:17.8 | r7_* -> 14.5:0.5:17.5
% If deriv is omitted, spline spec is reused for derivatives.

apply_s8_analysis_profile(analysis_profile, results_run);
slice_specs = s8_init_slice_specs(analysis_profile);

paths = get_project_paths();
if ~exist(paths.results_dir, 'dir')
    alt = s8_latest_results_run(paths.results_root);
    if ~isempty(alt)
        warning(['Results dir not found: %s\nUsing latest: %s'], paths.results_dir, ...
            fullfile(paths.results_root, alt));
        results_run = alt;
        apply_s8_analysis_profile(analysis_profile, results_run);
        paths = get_project_paths();
    end
end
pheno_xname = paths.pheno_xname;
pheno_linear = paths.pheno_linear_term;
if isempty(effect_vars)
    effect_vars = {pheno_xname, 'age'};
else
    % Remap stale pheno aliases so {'bmi','age'} still works under perc_median profiles.
    for i = 1:numel(effect_vars)
        if ismember(effect_vars{i}, {'bmi', 'perc_median'})
            effect_vars{i} = pheno_xname;
        end
    end
    effect_vars = unique(effect_vars, 'stable');
end
fprintf('s8 extended: profile=%s | pheno=%s | effect_vars=%s\n', ...
    analysis_profile, pheno_xname, strjoin(effect_vars, ', '));
fprintf('Results dir: %s\n', paths.results_dir);
fprintf('Shared BF dir: %s\n', paths.shared_datafiles);
fprintf('Pheno BF dir:  %s\n', paths.datafiles);

%% Design matrix and modality settings
% Combined-sex (r7_*, bmip_no_VE11): '1'. Sex-split full-sample/scanner: '1f' / '1m'.
if s8_combined_sex_profile(analysis_profile)
    desmatnum = '1';
else
    desmatnum = '1f';
end
modality = 'rni';
datatype = 'voxel';

dirname_out = paths.results_dir;
if ~exist(dirname_out, 'dir')
    error('Results directory not found: %s\nSet results_run at top of script or BMI_IMAGING_RESULTS_RUN in env.', dirname_out);
end

designmat_file = list_s8_designmat_dirs(dirname_out);
if isempty(designmat_file)
    error('No designmat* folders under %s', dirname_out);
end
fprintf('Designmat folders found: %s\n', strjoin(designmat_file, ', '));
dirname_results = fullfile(dirname_out, designmat_file);

addpath(genpath(paths.cmig_tools));
addpath(genpath(paths.matlab_utils));

%% Atlas
atlasVersion = 'ABCD3';
release = '6.0';
dirname_atlas = sprintf('/home/dpecheva/tools/internal/showVolData/Atlas/%s_%s_cor10', release, atlasVersion);
fname_cFA = sprintf('%s/CO_%s_%s_cor10.mat', dirname_atlas, release, atlasVersion);
vname_cFA = sprintf('CO_60_%s_cor10', atlasVersion);
load(fname_cFA);
fname_t1 = sprintf('%s/T1_%s_%s_cor10.mat', dirname_atlas, release, atlasVersion);
load(fname_t1);
vname_t1 = sprintf('T1_60_%s_cor10', atlasVersion);

%% Display limits (tune per modality / pheno)
pwaldlims_dict = struct(pheno_xname, [0, 30], 'age', [0, 100], 'pds', [0, 5]);
pwaldlims_W_dict = struct(pheno_xname, [0, 100], 'age', [0, 300], 'pds', [0, 50]);
pwaldlims_Bon_dict = pwaldlims_dict;
pwaldlims_FDR_dict = pwaldlims_dict;

%% Grayscale limits (showVol gray maps; empty struct / [] = use color limit above)
gwaldlims_dict = struct('bmi', [0 30], 'age', [0 100], 'pds', [0 5]);
gwaldlims_Bon_dict = struct('bmi', [0 5], 'age', [0 10], 'pds', [0 5]);
gwaldlims_FDR_dict = struct('bmi', [0 5], 'age', [0 10], 'pds', [0 5]);
gwaldlims_W_dict = struct('bmi', [0 100], 'age', [0 300], 'pds', [0 5]);

phenolim = [-0.0008 0.0008];
agelim = [-0.01 0.01];
pdslim = [-0.008 0.008];
sex_tlim = [-6 6];              % sex t = beta/SE display limits
derivlim_pheno = phenolim;
derivlim_age = agelim;
derivlim_pds = pdslim;

phenolim_gray = [-0.0008 0.0008];
agelim_gray = [-0.01 0.01];
pdslim_gray = [];
sex_tlim_gray = [];
derivlim_pheno_gray = [-0.0008 0.0008];
derivlim_age_gray = [-0.01 0.01];
derivlim_pds_gray = [];

showPval = false;
cmap_linear = redblackblue;
cmap_splines = redblackblue;
cmap_wald = cmap_linear(200:end, :);
bg_vol = [];
CSFprob = [];
interpMethod = [];

%% Basis functions and eval grids
dirname_basis = paths.datafiles;
shared_basis = paths.shared_datafiles;

tbl_bf_pheno = readtable(fullfile(dirname_basis, paths.pheno_bf_table));
tbl_bf_age = readtable(fullfile(shared_basis, 'table_bf_age.txt'));
tbl_bf_pds = readtable(fullfile(shared_basis, 'table_bf_pds.txt'));

vec_file = fullfile(dirname_basis, paths.pheno_vec_file);
readtable(vec_file, 'Delimiter', '\t'); %#ok<NASGU>
phenovec = pheno_bf_xvec(tbl_bf_pheno, paths);
if height(tbl_bf_pheno) ~= numel(phenovec)
    error('BF table rows (%d) != phenovec length (%d). Re-run s3.', height(tbl_bf_pheno), numel(phenovec));
end
fprintf('Pheno eval grid: n=%d, [%g : %g : %g]\n', numel(phenovec), phenovec(1), phenovec(2)-phenovec(1), phenovec(end));
agevec = pheno_bf_xvec(tbl_bf_age, paths, 'age');
fprintf('Age eval grid: n=%d, [%g : %g : %g]\n', numel(agevec), agevec(1), agevec(2)-agevec(1), agevec(end));
pdsvec = linspace(1, 5, 101);

fprintf('Resolving slice indices from slice_specs:\n');
pheno_spline_slices = resolve_s8_slices(phenovec, slice_specs, pheno_xname, 'spline');
pheno_deriv_slices = resolve_s8_slices(phenovec, slice_specs, pheno_xname, 'deriv');
age_spline_slices = resolve_s8_slices(agevec, slice_specs, 'age', 'spline');
age_deriv_slices = resolve_s8_slices(agevec, slice_specs, 'age', 'deriv');
pds_spline_slices = resolve_s8_slices(pdsvec, slice_specs, 'pds', 'spline');
pds_deriv_slices = resolve_s8_slices(pdsvec, slice_specs, 'pds', 'deriv');

abcd_sync_path = fema_default_abcd_sync_path();
dirname_imaging = fema_imaging_dirname(datatype, modality, abcd_sync_path);

des_idx = find(contains(designmat_file, sprintf('designmat%s_', desmatnum)))';
if isempty(des_idx)
    fprintf('Available designmat folders:\n');
    fprintf('  %s\n', designmat_file{:});
    error('No folder matches desmatnum=%s under %s', desmatnum, dirname_out);
end

for des = des_idx

    est_file = fema_results_path(dirname_results{des}, datatype, modality);
    R = fema_load_results('ResultsFile', est_file, ...
        'Datatype', datatype, 'FstemImaging', modality, 'DirnameImaging', dirname_imaging);

    beta_hat = R.beta_hat;
    coeffCovar = R.coeffCovar;
    mask = R.mask;
    colnames_model = R.colnames_model;
    if isfield(R, 'zmat')
        zmat = R.zmat;
    else
        zmat = [];
    end

    fprintf('s8 extended: loaded %s (datatype=%s modality=%s)\n', est_file, R.datatype, R.modality);

    Spheno_idx = find(find_cell(regexp(colnames_model, ['^', paths.pheno_bf_grep], 'match')));
    Sage_idx = find(find_cell(regexp(colnames_model, '^bf_age', 'match')));
    Spds_idx = find(find_cell(regexp(colnames_model, '^bf_pds', 'match')));
    sex_idx = find(contains(colnames_model, 'sexF'));
    intercept_idx = find(strcmp(colnames_model, 'intercept'));

    if ~strcmp(datatype, 'voxel')
        warning('s8 extended: skipping non-voxel datatype for %s', designmat_file{des});
        continue;
    end

    vols = [];

    % %% Sex linear t-stat from FEMA zmat (beta/SE), mapped with fullvol
    % if include_sex_tstat
    %     if isempty(sex_idx)
    %         fprintf('  No sexF column in %s; skipping sex t-stat map.\n', designmat_file{des});
    %     elseif isempty(zmat)
    %         warning('No zmat in %s; skipping sex t-stat map.', est_file);
    %     else
    %         vol_sex_t = fullvol(zmat(sex_idx(1), :), mask);
    %         statname = sprintf('t (sex linear), %s', designmat_file{des});
    %         vols = convert_and_append_vols(vols, vol_sex_t, modality, statname, ...
    %             colnames_model, sex_idx(1), sex_tlim, showPval, cmap_linear, bg_vol, ...
    %             CSFprob, interpMethod, atlasVersion, vol_display_mode, ...
    %             bmi_showvol_gray_limits(sex_tlim_gray, sex_tlim));
    %     end
    % end

    %% Wald chi-square and -log10(p) maps
    for i = 1:numel(effect_vars)
        xname = effect_vars{i};
        if ~isfield(pwaldlims_dict, xname)
            warning('No pwaldlims for xname=%s; skipping Wald maps.', xname);
            continue;
        end
        pwaldlimits = pwaldlims_dict.(xname);
        pwaldlimits_Bon = pwaldlims_Bon_dict.(xname);
        pwaldlimits_FDR = pwaldlims_FDR_dict.(xname);
        if isfield(pwaldlims_W_dict, xname)
            pwaldlimits_W = pwaldlims_W_dict.(xname);
        else
            pwaldlimits_W = pwaldlimits;
        end

        corrections = {};
        if include_logp
            corrections = wald_corrections;
        end

        vols = process_wald_stats(vols, dirname_results{des}, modality, xname, designmat_file{des}, ...
            pwaldlimits, pwaldlimits_Bon, pwaldlimits_FDR, showPval, cmap_wald, bg_vol, CSFprob, interpMethod, atlasVersion, ...
            'Corrections', corrections, ...
            'IncludeWaldStat', include_wald_stat, ...
            'PwaldlimitsW', pwaldlimits_W, ...
            'VolDisplayMode', vol_display_mode, ...
            'Graylimits', bmi_showvol_gray_limits(gwaldlims_dict, pwaldlimits, xname), ...
            'GraylimitsBon', bmi_showvol_gray_limits(gwaldlims_Bon_dict, pwaldlimits_Bon, xname), ...
            'GraylimitsFDR', bmi_showvol_gray_limits(gwaldlims_FDR_dict, pwaldlimits_FDR, xname), ...
            'PwaldlimitsWGray', bmi_showvol_gray_limits(gwaldlims_W_dict, pwaldlimits_W, xname));
    end

    %% Spline beta_hat slices and derivatives
    if any(contains(colnames_model, 'bf_'))
        spline_meta = {
            pheno_xname, tbl_bf_pheno, phenovec, Spheno_idx, phenolim, pheno_spline_slices, pheno_deriv_slices, derivlim_pheno;
            'age',       tbl_bf_age,   agevec,   Sage_idx,   agelim,   age_spline_slices,   age_deriv_slices,   derivlim_age;
            'pds',       tbl_bf_pds,   pdsvec,   Spds_idx,   pdslim,   pds_spline_slices,   pds_deriv_slices,   derivlim_pds;
        };

        for i = 1:size(spline_meta, 1)
            xname = spline_meta{i, 1};
            if ~ismember(xname, effect_vars)
                continue;
            end
            tbl_bf = spline_meta{i, 2};
            xvec = spline_meta{i, 3};
            idx_bf = spline_meta{i, 4};
            limits = spline_meta{i, 5};
            spline_slice_indices = spline_meta{i, 6};
            deriv_slice_indices = spline_meta{i, 7};
            deriv_limits = spline_meta{i, 8};

            if isempty(idx_bf)
                fprintf('  No spline columns for %s in %s; skipping splines/derivatives.\n', xname, designmat_file{des});
                continue;
            end

            if include_spline_slices
                vols = process_spline_slices(vols, tbl_bf, xvec, beta_hat, idx_bf, ...
                    intercept_idx, mask, coeffCovar, spline_slice_indices, ...
                    modality, designmat_file{des}, colnames_model, limits, ...
                    showPval, cmap_splines, bg_vol, CSFprob, interpMethod, atlasVersion, xname, vol_display_mode);
            end

            if include_derivatives
                vols = process_spline_derivatives(vols, tbl_bf, xvec, beta_hat, idx_bf, ...
                    intercept_idx, mask, coeffCovar, deriv_slice_indices, ...
                    modality, designmat_file{des}, colnames_model, deriv_limits, ...
                    showPval, cmap_splines, bg_vol, CSFprob, interpMethod, atlasVersion, xname, vol_display_mode);
            end
        end
    end

    coords = [0 -17 2];
    if isempty(vols)
        warning('s8 extended: no volumes to display for %s (check s7 Wald outputs).', designmat_file{des});
    else
        fprintf('s8 extended: opening showVol with %d volumes for %s\n', numel(vols), designmat_file{des});
        showVol(vols, eval(vname_cFA), eval(vname_t1), coords);
    end
end

function idx = resolve_s8_slices(xvec, slice_specs, xname, map_type)
%RESOLVE_S8_SLICES Pick spline or derivative slice indices for one predictor.
    spec = struct('mode', 'default', 'xname', xname);
    if isfield(slice_specs, xname)
        entry = slice_specs.(xname);
        if isfield(entry, map_type) && ~isempty(entry.(map_type))
            spec = entry.(map_type);
        elseif strcmp(map_type, 'deriv') && isfield(entry, 'spline') && ~isempty(entry.spline)
            spec = entry.spline;
        end
    end
    if ~isfield(spec, 'xname')
        spec.xname = xname;
    end
    idx = bmi_resolve_slice_indices(xvec, spec, sprintf('%s %s', xname, map_type));
end
