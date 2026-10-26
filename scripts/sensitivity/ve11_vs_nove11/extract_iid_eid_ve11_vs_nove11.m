%% Extract post-FEMA-intersect analytic iid/eid for VE11 vs noVE11.
%
% FEMA saves *_designMatrix.mat with the pre-intersect designmat rows (7,931/arm).
% The fit uses FEMA_intersect_design (designmat ∩ imaging QC: in vol_info and
% mean corrmat >= 0.8), which yields 7,889 VE11 and 7,579 noVE11.
%
% Writes iid_eid_rni.txt with the *post-intersect* visit list under each
% arm × pheno results designmat folder.
%
% Usage:
%   matlab -batch "run('scripts/matlab_add_paths.m'); run('scripts/sensitivity/ve11_vs_nove11/extract_iid_eid_ve11_vs_nove11.m')"
%
% Defaults:
%   bmip        -> results_2026_09_02 / designmats_2026_09_02
%   perc_median -> results_2026_10_02 / designmats_2026_10_02

clearvars -except pheno_list arm_list results_run_by_pheno modality ...
    corrvec_thresh dirname_imaging;
close all;

%% Knobs
pheno_list = {'bmip', 'perc_median'};
arm_list = {'VE11', 'noVE11'};
modality = 'rni';
corrvec_thresh = 0.8;
results_run_by_pheno = struct( ...
    'bmip', 'results_2026_09_02', ...
    'perc_median', 'results_2026_10_02');

run(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'matlab_add_paths.m'));

if ~exist('dirname_imaging', 'var') || isempty(dirname_imaging)
    sync = strtrim(getenv('BMI_IMAGING_ABCD_SYNC'));
    if isempty(sync)
        sync = '/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0';
    end
    dirname_imaging = fullfile(sync, 'concatenated', 'imaging', 'voxelwise', 'rsi');
end

[pid_vi, sid_vi, corrvec] = local_load_volinfo_corr(dirname_imaging);
pass = corrvec >= corrvec_thresh;
keys_ok = string(strcat(string(pid_vi(pass)), '_', string(sid_vi(pass))));
fprintf('Imaging QC pass (corr>=%.2f): %d visits under %s\n', ...
    corrvec_thresh, numel(keys_ok), dirname_imaging);

for ip = 1:numel(pheno_list)
    pheno = pheno_list{ip};
    results_run = results_run_by_pheno.(pheno);
    designmats_tag = strrep(results_run, 'results_', 'designmats_');
    for ia = 1:numel(arm_list)
        arm = arm_list{ia};
        profile = sprintf('ve11_nove11_%s_%s', arm, pheno);
        apply_s8_analysis_profile(profile, results_run);
        paths = get_project_paths();
        results_dir = paths.results_dir;
        arm_root = fileparts(fileparts(results_dir)); % .../<pheno>/<arm>

        d = dir(fullfile(results_dir, 'designmat1_*'));
        models = {d([d.isdir]).name};
        models = models(startsWith(models, 'designmat1_'));
        if isempty(models)
            warning('No designmat1_* under %s', results_dir);
            continue;
        end

        for j = 1:numel(models)
            model = models{j};
            dir_results = fullfile(results_dir, model);
            des_txt = fullfile(arm_root, 'designmats', designmats_tag, [model, '.txt']);
            if ~isfile(des_txt)
                warning('Designmat not found: %s', des_txt);
                continue;
            end

            T = readtable(des_txt);
            [keys, pid, sid] = local_visit_keys(T);
            keep = ismember(string(keys), keys_ok);

            iid = string(pid(keep));
            eid = string(sid(keep));
            t = table(iid(:), eid(:), 'VariableNames', {'iid', 'eid'});
            output_file = fullfile(dir_results, sprintf('iid_eid_%s.txt', modality));
            writetable(t, output_file);
            fprintf('\n=== %s / %s / %s ===\n', pheno, arm, model);
            fprintf('  designmat=%d  post-intersect=%d  dropped=%d\n', ...
                numel(keys), height(t), numel(keys) - height(t));
            fprintf('  Wrote %s\n', output_file);
        end
    end
end

fprintf('\nextract_iid_eid_ve11_vs_nove11: done (post-FEMA-intersect lists).\n');

function [pid, sid, corrvec] = local_load_volinfo_corr(dirname_imaging)
fname = fullfile(dirname_imaging, 'vol_info.mat');
S = load(fname, 'participant_id', 'session_id', 'corrmat');
pid = cellstr(string(S.participant_id(:)));
sid = cellstr(string(S.session_id(:)));
cm = S.corrmat;
if size(cm, 1) < size(cm, 2)
    corrvec = mean(cm, 1);
else
    corrvec = mean(cm, 2);
end
corrvec = corrvec(:);
end

function [keys, pid, sid] = local_visit_keys(T)
if all(ismember({'participant_id', 'session_id'}, T.Properties.VariableNames))
    pid = cellstr(string(T.participant_id));
    sid = cellstr(string(T.session_id));
elseif all(ismember({'iid', 'eid'}, T.Properties.VariableNames))
    pid = cellstr(string(T.iid));
    sid = cellstr(string(T.eid));
else
    error('Table needs participant_id/session_id or iid/eid.');
end
keys = cellstr(strcat(string(pid), '_', string(sid)));
end
