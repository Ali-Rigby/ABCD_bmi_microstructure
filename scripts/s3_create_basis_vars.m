%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basis functions: pheno (bmip / perc_median) + shared age + PDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

paths = get_project_paths();
addpath(genpath(paths.cmig_tools));
addpath(paths.matlab_utils);

pheno = getenv('BMI_PHENO_VARIABLE');
switch pheno
    case 'bmip'
        phenoLabel = 'BMI';
        vecFile = 'bmi_percentile_vector.txt';
        outRaw = 'table_bf_bmi_perc.txt';
        outInterp = 'table_bf_bmi_perc_interp.txt';
        phenoKnots = @(v) prctile(v, [25 50 75]);
    case 'perc_median'
        phenoLabel = 'PercMedian';
        vecFile = 'perc_median_vector.txt';
        outRaw = 'table_bf_perc_median.txt';
        outInterp = 'table_bf_perc_median_interp.txt';
        phenoKnots = @(v) prctile(v, [25 50 75]);
    otherwise
        error('Set BMI_PHENO_VARIABLE to bmip or perc_median.');
end

cfg = struct();
cfg.splineType = 'nsk';
cfg.Xpowers = [0];
cfg.method = 'svd';
cfg.outDir = [];
cfg.intercept = false;
cfg.optCommand = '';
cfg.optAppend = '';
cfg.cleanUp = true;
cfg.instance = 1;

fprintf('\n--- Basis functions (unified BMI imaging) ---\n');
fprintf('Pheno: %s | profile datafiles: %s | shared: %s\n', pheno, paths.datafiles, paths.shared_datafiles);

cfg.basepath = paths.datafiles;
val_table = readtable(fullfile(cfg.basepath, vecFile), 'Delimiter', '\t');
valvec = val_table.Var1;
phenoMinMax = pheno_basis_minmax(valvec, paths);
if strcmp(pheno, 'bmip')
    fprintf('Pheno eval grid (bmip fixed): [%g, %g], n=%d\n', ...
        phenoMinMax(1), phenoMinMax(2), phenoMinMax(3));
else
    fprintf('Pheno eval grid (BMI_IMAGING_PHENO_GRID_STEP): [%g, %g], n=%d\n', ...
        phenoMinMax(1), phenoMinMax(2), phenoMinMax(3));
end

run_basis_job(cfg, phenoLabel, vecFile, phenoMinMax, outRaw, outInterp, phenoKnots, @createBasisFunctions);

% Shared age/PDS (shared/datafiles; rebuilt every s3 run)
sharedCfg = cfg;
sharedCfg.basepath = paths.shared_datafiles;
run_basis_job(sharedCfg, 'Age', 'age_vector_from_data.txt', [], 'table_bf_age.txt', 'table_bf_age_interp.txt', ...
    @(~) [9.96 12.0 14.1 16.1], @createBasisFunctions_round);
run_basis_job(sharedCfg, 'PDS', 'pds_vector_from_data.txt', [1 5 101], 'table_bf_pds.txt', 'table_bf_pds_interp.txt', ...
    @(~) [1 2 3 4 5], @createBasisFunctions);

fprintf('\n--- Done ---\n');

function run_basis_job(cfg, label, input_file, minMax, output_raw, output_interp, knots_fn, basis_fn)
fprintf('\n--- Processing %s ---\n', label);
val_table = readtable(fullfile(cfg.basepath, input_file), 'Delimiter', '\t');
valvec = val_table.Var1;
fprintf('Length of loaded %s vector: %d\n', label, length(valvec));
knots = knots_fn(valvec);
fprintf('Knots: %s\n', mat2str(knots));
if isempty(minMax)
    minX = round(min(valvec) * 10) / 10;
    maxX = ceil(max(valvec) * 10) / 10;
    minMax = [minX maxX 101];
end
[bfs_interp, bfs_raw, ~, ~, bf_settings, ~] = basis_fn( ...
    valvec, knots, cfg.splineType, cfg.Xpowers, cfg.method, cfg.outDir, minMax, cfg.intercept, ...
    cfg.optCommand, cfg.optAppend, cfg.cleanUp, cfg.instance);
if isfield(bf_settings, 'knots')
    fprintf('Final knots: %s\n', mat2str(bf_settings.knots(:)', 4));
end
fprintf('Rows in %s basis matrix: %d\n', label, size(bfs_interp, 1));
writematrix(bfs_raw, fullfile(cfg.basepath, output_raw), 'Delimiter', '\t');
writematrix(bfs_interp, fullfile(cfg.basepath, output_interp), 'Delimiter', '\t');
end
