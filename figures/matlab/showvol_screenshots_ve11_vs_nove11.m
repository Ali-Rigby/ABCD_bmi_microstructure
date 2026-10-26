%% Save showVol screenshots for VE11 vs noVE11 (matched 6.0) bundles.
% Sibling of showVol_screenshots.m. Prerequisite:
%   run('figures/matlab/showvol_results_for_screenshots_ve11_vs_nove11.m');
%
% Loops pheno × arm × summary modality (defaults cover all 4 sample arms).
%
% GUI: edit knobs below, then run.
% Shell:
%   matlab -batch "run('scripts/matlab_addpaths.m'); run('figures/matlab/showvol_screenshots_ve11_vs_nove11.m')"

clearvars -except pheno_list arm_list desmat_filter summary_modalities ...
    results_run_by_pheno;
close all;

%% User knobs
desmat_filter = '1';   % combined-sex designmat1_*
pheno_list    = {'bmip', 'perc_median'};
arm_list      = {'VE11', 'noVE11'};
% One or both bundle types written by the prep sibling script
summary_modalities = {'rni_wald', 'rni_derivatives'};

results_run_by_pheno = struct( ...
    'bmip', 'results_2026_09_02', ...
    'perc_median', 'results_2026_10_02');

%% Paths
run(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'scripts', 'matlab_add_paths.m'));
addpath(genpath('/usr/pubsw/packages/cmig/cmig_tools'));

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

coords = [0 -17 2];

for ip = 1:numel(pheno_list)
    pheno = strtrim(pheno_list{ip});
    if ~isfield(results_run_by_pheno, pheno)
        error('Add results_run_by_pheno.%s', pheno);
    end
    results_run = results_run_by_pheno.(pheno);

    for ia = 1:numel(arm_list)
        arm = strtrim(arm_list{ia});
        profile = sprintf('ve11_nove11_%s_%s', arm, pheno);
        apply_s8_analysis_profile(profile, results_run);

        paths = get_project_paths();
        datatype = 'voxel';
        xname = paths.pheno_xname;

        dirname_base = paths.results_dir;
        dirname_input_root = fullfile(dirname_base, ...
            sprintf('results_for_screenshots_%s_%s', datatype, xname));

        if ~exist(dirname_input_root, 'dir')
            warning('Missing screenshot bundle root: %s\n  Run showVol_results_for_screenshots_ve11_vs_nove11.m first.', ...
                dirname_input_root);
            continue;
        end

        designmat_file = list_s8_designmat_dirs(dirname_input_root);
        if isempty(designmat_file)
            warning('No designmat* under %s', dirname_input_root);
            continue;
        end
        dirname_input_results = fullfile(dirname_input_root, designmat_file);
        des_idx = select_designmat_indices(designmat_file, desmat_filter);
        if isempty(des_idx)
            warning('No designmat folders match filter ''%s'' under %s', ...
                desmat_filter, dirname_input_root);
            continue;
        end

        for im = 1:numel(summary_modalities)
            summary_modality = strtrim(summary_modalities{im});

            fprintf('\n=== showVol screenshots: %s / %s | %s ===\n', ...
                pheno, arm, summary_modality);
            fprintf('Results dir: %s\n', dirname_base);
            fprintf('Input root: %s\n', dirname_input_root);

            for des = des_idx(:)'
                results_path = dirname_input_results{des};
                results_file = dir(fullfile(results_path, ...
                    sprintf('vols_summary_%s.mat', summary_modality)));

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

                hfig = showVol(vols, coords, 1);
                drawnow;
                pause(3);

                fstem = labels;

                CoordSpace = 'RCS';
                coords = [97 101 146];

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

                save_showVol_images_ali(coords, outdir, fstem, ...
                    'drawROI', drawROI, 'parcellation', parcellation, ...
                    'roiNames', roiNames, 'CoordSpace', CoordSpace, ...
                    'lw', lw, 'fixedColor', fixedColor, ...
                    'hideColorbar', true, ...
                    'DesignmatName', DesignmatName, 'Modality', Modality);

                close(hfig);

                anatomical_mods = {vname_t1, vname_cFA};
                anatomical_labels = {'T1', 'cFA'};

                for anat_idx = 1:numel(anatomical_mods)
                    vname = anatomical_mods{anat_idx};
                    anat_label = anatomical_labels{anat_idx};

                    if exist(vname, 'var')
                        anat_vol = eval(vname);
                        hfig_anat = showVol({anat_vol}, coords, 1);
                        drawnow;
                        pause(3);

                        fstem_anat = {sprintf('%s_%s_%s_%s', anat_label, pheno, arm, desmat_tag)};
                        save_showVol_images_ali(coords, outdir, fstem_anat, ...
                            'drawROI', drawROI, 'parcellation', parcellation, ...
                            'roiNames', roiNames, 'CoordSpace', CoordSpace, ...
                            'lw', lw, 'fixedColor', fixedColor, ...
                            'hideColorbar', true, ...
                            'DesignmatName', DesignmatName, 'Modality', anat_label);

                        close(hfig_anat);
                    else
                        warning('Anatomical volume %s not found. Skipping.', vname);
                    end
                end
            end
        end
    end
end

fprintf('\nDone: VE11 vs noVE11 showVol screenshots.\n');

%% Local helpers
function idx = select_designmat_indices(designmat_names, desmat_filter)
    if strcmpi(strtrim(desmat_filter), 'all')
        idx = (1:numel(designmat_names))';
        return;
    end
    pattern = sprintf('designmat%s_', desmat_filter);
    idx = find(contains(designmat_names, pattern))';
end

function tag = designmat_tag_from_name(dname)
    tok = regexp(dname, '^designmat(\d+[fm]?)_', 'tokens', 'once');
    if isempty(tok)
        tag = 'unknown';
    else
        tag = tok{1};
    end
end
