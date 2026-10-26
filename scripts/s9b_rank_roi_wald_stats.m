% Rank ROIs by mean Wald statistic from s7 voxelwise outputs.
%
% Reads s7 <modality>_Wald_splines_CI_output_<xname>.mat (vol_W_s by default),
% averages uncorrected Wald stats across all in-mask voxels in each ROI
% (not FDR/Bonferroni filtered), counts FDR-significant voxels per ROI
% (vol_pW_*_FDR == 0), and writes ranked CSVs.
%
% ROI groups (see utilities/matlab_utils/bmi_roi_wald_config.m):
%   bilateral - aseg left+right combined
%   fiber       - Fmin, Fmaj, Bilateral-CST
%   pauli       - Pauli subcortical nuclei
%
% Before running (MATLAB GUI or batch):
%   setenv('BMI_PHENO_VARIABLE', 'bmip');          % or perc_median
%   setenv('BMI_IMAGING_ANALYSIS_DIR', 'bmip');
%   setenv('BMI_IMAGING_CEPH_ROOT', '/space/ceph/1/ABCD/users/arigby/2026_bmi_imaging_analysis');
%   setenv('BMI_IMAGING_RESULTS_RUN', 'results_2026_06_29');
%   setenv('BMI_IMAGING_PHENO_GRID_STEP', '1');    % required for perc_median
%
%   run('scripts/matlab_add_paths.m');
%   s9b_rank_roi_wald_stats
%
% Optional env:
%   BMI_IMAGING_S7_FSTEM=rni
%   BMI_ROI_WALD_FIELD=vol_W_s   % vol_W_l | vol_W_s | vol_W_ls
%   BMI_ROI_WALD_MIN_VOXELS=10

paths = get_project_paths();
addpath(genpath(paths.cmig_tools));
addpath(paths.matlab_utils);

fprintf('s9b rank ROI Wald: BMI_PHENO_VARIABLE=%s\n', paths.pheno_variable);
fprintf('  results_dir = %s\n', paths.results_dir);

modality = getenv('BMI_IMAGING_S7_FSTEM');
if isempty(modality), modality = 'rni'; end

wald_field = getenv('BMI_ROI_WALD_FIELD');
if isempty(wald_field), wald_field = 'vol_W_s'; end

% Map Wald field -> matching FDR indicator volume from s7
switch wald_field
    case 'vol_W_s'
        fdr_field = 'vol_pW_s_FDR';
    case 'vol_W_l'
        fdr_field = 'vol_pW_l_FDR';
    case 'vol_W_ls'
        fdr_field = 'vol_pW_ls_FDR';
    otherwise
        fdr_field = '';
        warning('No FDR field mapping for %s; n_fdr_sig will be 0.', wald_field);
end

min_voxels = str2double(getenv('BMI_ROI_WALD_MIN_VOXELS'));
if isnan(min_voxels) || min_voxels < 1
    min_voxels = 1;
end

xname = paths.pheno_xname;
cfg = bmi_roi_wald_config();
datatype = 'voxel';
abcd_sync_path = fema_default_abcd_sync_path();

out_root = fullfile(paths.results_dir, 'roi_wald_rankings');
if ~exist(out_root, 'dir')
    mkdir(out_root);
end

designmat_dirs = dir(paths.results_dir);
designmat_dirs = designmat_dirs([designmat_dirs(:).isdir]);
names = {designmat_dirs(:).name};
designmat_dirs = designmat_dirs(startsWith(names, 'designmat'));

all_tables = table();

for d = 1:numel(designmat_dirs)
    dirname_desmat = designmat_dirs(d).name;
    dirname_out = fullfile(paths.results_dir, dirname_desmat);
    fprintf('\n--- %s ---\n', dirname_desmat);

  if contains(dirname_desmat, 'female')
    sex_label = 'female';
  elseif contains(dirname_desmat, 'male')
    sex_label = 'male';
  else
    sex_label = 'unknown';
  end

    wald_file = fullfile(dirname_out, sprintf('%s_Wald_splines_CI_output_%s.mat', modality, xname));
    if ~isfile(wald_file)
        warning('Wald file not found: %s', wald_file);
        continue;
    end

    S = load(wald_file);
    if ~isfield(S, wald_field)
        warning('Field %s missing in %s', wald_field, wald_file);
        continue;
    end
    vol_W = S.(wald_field);

    vol_FDR = [];
    if ~isempty(fdr_field)
        if isfield(S, fdr_field)
            vol_FDR = S.(fdr_field);
        else
            warning('FDR field %s missing in %s; n_fdr_sig will be 0.', fdr_field, wald_file);
        end
    end

    est_file = fema_results_path(dirname_out, datatype, modality);
    if ~isfile(est_file)
        warning('FEMA results not found: %s', est_file);
        continue;
    end
    dirname_imaging = fema_imaging_dirname(datatype, modality, abcd_sync_path);
    R = fema_load_results('ResultsFile', est_file, ...
        'Datatype', datatype, 'FstemImaging', modality, 'DirnameImaging', dirname_imaging);
    brain_mask = R.mask;

    T = bmi_rank_roi_wald(vol_W, brain_mask, cfg, min_voxels, vol_FDR);
    if isempty(T)
        warning('No ROIs ranked for %s', dirname_desmat);
        continue;
    end

    T.sex = repmat({sex_label}, height(T), 1);
    T.designmat = repmat({dirname_desmat}, height(T), 1);
    T.pheno = repmat({paths.pheno_variable}, height(T), 1);
    T.xname = repmat({xname}, height(T), 1);
    T.wald_field = repmat({wald_field}, height(T), 1);

    disp(T);

    out_csv = fullfile(out_root, sprintf('%s_%s_%s_%s_ranked_wald.csv', modality, xname, sex_label, wald_field));
    writetable(T, out_csv);
    fprintf('Saved: %s\n', out_csv);

    all_tables = [all_tables; T]; %#ok<AGROW>
end

if ~isempty(all_tables)
    combined_csv = fullfile(out_root, sprintf('%s_%s_%s_all_sexes_ranked_wald.csv', modality, xname, wald_field));
    writetable(all_tables, combined_csv);
    fprintf('\nCombined rankings: %s\n', combined_csv);
else
    fprintf('\nNo rankings produced.\n');
end

fprintf('\nDone. Outputs under: %s\n', out_root);
