%% Save showVol screenshots from pre-built volume summary bundles.
% Run showvol_results_for_screenshots.m first.
%
% GUI-friendly: set pheno_list and results_run below (no .csh required).
% Atlas path: BMI_IMAGING_SHOWVOL_ATLAS or paths.showvol_atlas from get_project_paths.

%% User knobs
desmat_filter = '1f';              % '1f' | '1m' | '1' | 'all'
summary_modality = 'rni_derivatives';  % 'rni_derivatives' or 'rni_wald'
pheno_list = {'bmip'};             % {} = use BMI_PHENO_VARIABLE from env
results_run = 'results_2026_06_29'; % '' = env / auto date (see get_project_paths)

%% Paths (this file lives in figures/matlab/)
this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
run(fullfile(repo_root, 'scripts', 'matlab_add_paths.m'));

if isempty(pheno_list)
    pheno_env = strtrim(getenv('BMI_PHENO_VARIABLE'));
    if isempty(pheno_env)
        error(['Set BMI_PHENO_VARIABLE or pheno_list in this script.']);
    end
    pheno_list = {pheno_env};
end

%% Atlas underlays (set once; reused across phenotypes)
% Bootstrap paths for cmig_tools + atlas location (pheno set temporarily if needed)
if isempty(getenv('BMI_PHENO_VARIABLE')) && ~isempty(pheno_list)
    setenv('BMI_PHENO_VARIABLE', pheno_list{1});
    setenv('BMI_IMAGING_ANALYSIS_DIR', pheno_list{1});
end
paths0 = get_project_paths();
addpath(genpath(paths0.cmig_tools));

atlasVersion = 'ABCD3';
release = '6.0';
dirname_atlas = paths0.showvol_atlas;
fname_cFA = fullfile(dirname_atlas, sprintf('CO_%s_%s_cor10.mat', release, atlasVersion));
vname_cFA = sprintf('CO_60_%s_cor10', atlasVersion);
load(fname_cFA);
fname_t1 = fullfile(dirname_atlas, sprintf('T1_%s_%s_cor10.mat', release, atlasVersion));
load(fname_t1);
vname_t1 = sprintf('T1_60_%s_cor10', atlasVersion);

coords = [0 -17 2];

for ip = 1:numel(pheno_list)
    pheno = strtrim(pheno_list{ip});
    setenv('BMI_PHENO_VARIABLE', pheno);
    setenv('BMI_IMAGING_ANALYSIS_DIR', pheno);
    if ~isempty(results_run)
        setenv('BMI_IMAGING_RESULTS_RUN', results_run);
    end

    paths = get_project_paths();
    datatype = 'voxel';
    xname = paths.pheno_xname;

    dirname_base = paths.results_dir;
    dirname_input_root = fullfile(dirname_base, ...
        sprintf('results_for_screenshots_%s_%s', datatype, xname));

    if ~exist(dirname_input_root, 'dir')
        warning('Missing screenshot bundle root: %s (run showvol_results_for_screenshots first)', ...
            dirname_input_root);
        continue;
    end

    designmat_file = list_designmat_dirs(dirname_input_root);
    dirname_input_results = fullfile(dirname_input_root, designmat_file);
    des_idx = select_designmat_indices(designmat_file, desmat_filter);

    fprintf('\n=== showVol screenshots: %s | %s ===\n', pheno, summary_modality);
    fprintf('Results dir: %s\n', dirname_base);
    fprintf('Input root: %s\n', dirname_input_root);

    for des = des_idx
        results_path = dirname_input_results{des};
        results_file = dir(fullfile(results_path, sprintf('vols_summary_%s.mat', summary_modality)));

        if isempty(results_file)
            warning('No summary file vols_summary_%s.mat in %s. Skipping.', ...
                summary_modality, results_path);
            continue;
        end

        fname_results = fullfile(results_path, results_file(1).name);
        S = load(fname_results);
        if ~isfield(S, 'vols') || ~isfield(S, 'labels')
            error('Expected vols and labels in %s', fname_results);
        end
        vols = S.vols;
        labels = S.labels;

        outdir = dirname_base;
        if ~exist(outdir, 'dir')
            mkdir(outdir);
        end

        hfig = showVol(vols, coords, 3);
        drawnow;
        pause(3);

        fstem = labels;
        Labels = {};

        CoordSpace = 'RCS';
        coords = [69 100 147];

        drawROI = true;
        parcellation = {'aseg', 'fiber'};
        roiNames = struct( ...
            'aseg', {{'Left-Thalamus-Proper', 'Right-Thalamus-Proper', ...
                    'Left-Caudate', 'Right-Caudate', ...
                    'Left-Putamen', 'Right-Putamen', ...
                    'Left-Pallidum', 'Right-Pallidum', ...
                    'Left-Accumbens-area', 'Right-Accumbens-area'}}, ...
            'fiber', {{'R_ATR', 'L_ATR', 'R_CgC', 'L_CgC', 'R_CgH', 'L_CgH', 'Fmin'}});

        lw = 2.5;
        fixedColor = [0, 0.7059, 0.3529];

        rawDesignmatFolder = designmat_file{des};
        DesignmatName = rawDesignmatFolder;
        Modality = summary_modality;
        desmat_tag = designmat_tag_from_name(rawDesignmatFolder);

        save_showvol_images(coords, outdir, fstem, ...
            'drawROI', drawROI, 'parcellation', parcellation, ...
            'roiNames', roiNames, 'CoordSpace', CoordSpace, ...
            'lw', lw, 'fixedColor', fixedColor, ...
            'hideColorbar', false, ...
            'DesignmatName', DesignmatName, 'Modality', Modality);

        close(hfig);

        anatomical_mods = {vname_t1, vname_cFA};
        anatomical_labels = {'T1', 'cFA'};

        for anat_idx = 1:numel(anatomical_mods)
            vname = anatomical_mods{anat_idx};
            anat_label = anatomical_labels{anat_idx};

            if exist(vname, 'var')
                anat_vol = eval(vname);
                hfig_anat = showVol({anat_vol}, coords, 3);
                drawnow;
                pause(3);

                fstem_anat = {sprintf('%s_%s', anat_label, desmat_tag)};
                save_showvol_images(coords, outdir, fstem_anat, ...
                    'drawROI', drawROI, 'parcellation', parcellation, ...
                    'roiNames', roiNames, 'CoordSpace', CoordSpace, ...
                    'lw', lw, 'fixedColor', fixedColor, ...
                    'hideColorbar', false, ...
                    'DesignmatName', DesignmatName, 'Modality', anat_label);

                close(hfig_anat);
            else
                warning('Anatomical volume %s not found. Skipping.', vname);
            end
        end
    end
end

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
