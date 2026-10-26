% Wald stats for unified BMI imaging pipeline (s(bmip) or s(perc_median)).

clear all; close all;

paths = get_project_paths();
fprintf('s7 (pheno): BMI_PHENO_VARIABLE=%s\n', paths.pheno_variable);

addpath(genpath(paths.cmig_tools))
addpath(genpath(paths.matlab_utils))

% Default imaging identity when results are bare FEMA_estimates.mat (no outPrefix).
default_datatype = getenv('BMI_IMAGING_S7_DATATYPE');
if isempty(default_datatype), default_datatype = 'voxel'; end
default_fstem = getenv('BMI_IMAGING_S7_FSTEM');
if isempty(default_fstem), default_fstem = 'rni'; end

modalities_filter = {'rni'};
env_modalities = getenv('BMI_IMAGING_S7_MODALITIES');
if ~isempty(env_modalities)
    modalities_filter = strsplit(env_modalities, ',');
    modalities_filter = cellfun(@strtrim, modalities_filter, 'UniformOutput', false);
    modalities_filter = modalities_filter(~cellfun(@isempty, modalities_filter));
end
if isempty(modalities_filter)
    fprintf('s7: modalities_filter is empty — processing all modalities.\n');
    modalities_filter_lower = {};
else
    fprintf('s7: restricting to modality tag(s): %s\n', strjoin(modalities_filter, ', '));
    modalities_filter_lower = cellfun(@lower, modalities_filter, 'UniformOutput', false);
end

dirname_outroot = paths.results_dir;
fprintf('s7: results_dir = %s\n', dirname_outroot);
fprintf('s7: datafiles_dir = %s\n', paths.datafiles);

designmat_dirs = dir(dirname_outroot);
designmat_dirs = designmat_dirs([designmat_dirs(:).isdir]);
names = {designmat_dirs(:).name};
designmat_dirs = designmat_dirs(startsWith(names, 'designmat'));

bf_dir = paths.datafiles;
bf_dir_shared = paths.shared_datafiles;
linear_terms = {paths.pheno_linear_term, 'ab_g_dyn__visit_age', 'pds_categ_y_p_all'};
xnames = {paths.pheno_xname, 'age', 'pds'};

% Optional: restrict Wald/spline terms, e.g. BMI_IMAGING_S7_XNAMES=bmi
% (comma-separated xnames). Useful for XA30 where age/PDS BFs were rank-dropped.
env_xnames = strtrim(getenv('BMI_IMAGING_S7_XNAMES'));
if ~isempty(env_xnames)
    keep = strsplit(env_xnames, ',');
    keep = cellfun(@strtrim, keep, 'UniformOutput', false);
    keep = keep(~cellfun(@isempty, keep));
    keep_mask = ismember(xnames, keep);
    if ~any(keep_mask)
        error('BMI_IMAGING_S7_XNAMES=%s matched none of: %s', env_xnames, strjoin(xnames, ', '));
    end
    linear_terms = linear_terms(keep_mask);
    xnames = xnames(keep_mask);
    fprintf('s7: restricting Wald terms to xnames: %s\n', strjoin(xnames, ', '));
end

for d = 1:length(designmat_dirs)
    dirname_desmat = designmat_dirs(d).name;
    dirname_out = fullfile(dirname_outroot, dirname_desmat);
    fprintf('Processing design matrix folder: %s\n', dirname_desmat);

    runs = fema_list_results(dirname_out, ...
        'DefaultDatatype', default_datatype, 'DefaultFstem', default_fstem);

    for r = 1:numel(runs)
        run = runs(r);
        modality = run.modality;

        if ~isempty(modalities_filter_lower)
            if ~ismember(lower(modality), modalities_filter_lower)
                fprintf('  Skipping %s (modality %s not in filter)\n', run.fname_results, modality);
                continue;
            end
        end

        fprintf('  Processing %s (datatype=%s modality=%s)\n', ...
            run.fname_results, run.datatype, modality);

        switch lower(run.datatype)
            case {'voxel', 'voxelwise'}
                fprintf('    -> Detected voxel-based data\n');
                process_voxel_data(run.fname_results, dirname_out, linear_terms, xnames, ...
                    dirname_desmat, bf_dir, ...
                    'Datatype', run.datatype, 'FstemImaging', run.fstem_imaging, ...
                    'BfDirShared', bf_dir_shared);
            case {'vertex', 'vertexwise'}
                fprintf('    -> Detected vertex-based data (processVertexData not implemented)\n');
            otherwise
                fprintf('    -> Unknown datatype ''%s'', skipping\n', run.datatype);
        end
    end
end
