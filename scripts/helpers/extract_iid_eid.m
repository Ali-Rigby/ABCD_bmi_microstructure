% extract_iid_eid.m
% Write *post-FEMA_intersect_design* visit lists (iid/eid) for downstream demo tables.
%
% FEMA 4.x *_designMatrix.mat keeps pre-intersect designmat rows. The fit uses
% designmat ∩ imaging QC (in vol_info and mean corrmat >= 0.8). This script
% rebuilds that intersect and writes iid_eid_<modality>.txt under each
% designmat1f_/designmat1m_ results folder.
%
% Usage (after s6):
%   source env/run_bmip.env        % or run_perc_median.env
%   export BMI_IMAGING_RESULTS_RUN=results_2026_06_29
%   export BMI_IMAGING_RUN_TAG=2026_06_29
%   export BMI_IMAGING_AUTO_RUN_DATE=0
%   matlab -batch "run('scripts/matlab_add_paths.m'); run('scripts/helpers/extract_iid_eid.m')"
%
% Env:
%   BMI_IMAGING_EXTRACT_MODALITIES  comma-separated fstem tags (default: rni)
%   BMI_IMAGING_S7_DATATYPE         unused here (kept for compatibility)
%   BMI_IMAGING_ABCD_SYNC           imaging root (default: parekh 6.0 concat)
%   BMI_IMAGING_CORRVEC_THRESH      default 0.8

paths = get_project_paths();
fprintf('extract_iid_eid: pheno=%s results_dir=%s\n', paths.pheno_variable, paths.results_dir);

modalities_env = strtrim(getenv('BMI_IMAGING_EXTRACT_MODALITIES'));
if isempty(modalities_env)
    modalities = {'rni'};
else
    modalities = strsplit(modalities_env, ',');
    modalities = cellfun(@strtrim, modalities, 'UniformOutput', false);
    modalities = modalities(~cellfun(@isempty, modalities));
end

corr_env = strtrim(getenv('BMI_IMAGING_CORRVEC_THRESH'));
if isempty(corr_env)
    corrvec_thresh = 0.8;
else
    corrvec_thresh = str2double(corr_env);
end

sync = strtrim(getenv('BMI_IMAGING_ABCD_SYNC'));
if isempty(sync)
    sync = '/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0';
end
dirname_imaging = fullfile(sync, 'concatenated', 'imaging', 'voxelwise', 'rsi');
vox_override = strtrim(getenv('BMI_IMAGING_VOXELWISE_RSI_DIR'));
if ~isempty(vox_override)
    dirname_imaging = vox_override;
end

results_dir = paths.results_dir;
if ~isfolder(results_dir)
    error('extract_iid_eid: results_dir not found: %s', results_dir);
end

results_run = paths.results_run;
designmats_tag = strrep(results_run, 'results_', 'designmats_');
designmats_dir = fullfile(paths.analysis_base, 'designmats', designmats_tag);
if ~isfolder(designmats_dir) && isfolder(paths.designmats_dir)
    designmats_dir = paths.designmats_dir;
end
if ~isfolder(designmats_dir)
    error('extract_iid_eid: designmats dir not found: %s', designmats_dir);
end

[pid_vi, sid_vi, corrvec] = local_load_volinfo_corr(dirname_imaging);
pass = corrvec >= corrvec_thresh;
keys_ok = string(strcat(string(pid_vi(pass)), '_', string(sid_vi(pass))));
fprintf('Imaging QC pass (corr>=%.2f): %d visits under %s\n', ...
    corrvec_thresh, numel(keys_ok), dirname_imaging);
fprintf('Designmats: %s\n', designmats_dir);

d = dir(fullfile(results_dir, 'designmat*'));
models = {d([d.isdir]).name};
models = models(startsWith(models, 'designmat1f_') | startsWith(models, 'designmat1m_'));

if isempty(models)
    error('extract_iid_eid: no designmat1f_/designmat1m_ folders under %s', results_dir);
end

fprintf('Found %d designmat result folder(s)\n', numel(models));

for j = 1:numel(models)
    model = models{j};
    dir_results = fullfile(results_dir, model);
    des_txt = fullfile(designmats_dir, [model, '.txt']);
    if ~isfile(des_txt)
        warning('extract_iid_eid: designmat not found: %s', des_txt);
        continue;
    end

    T = readtable(des_txt);
    [keys, pid, sid] = local_visit_keys(T);
    keep = ismember(string(keys), keys_ok);
    iid_keep = string(pid(keep));
    eid_keep = string(sid(keep));

    for i = 1:numel(modalities)
        modality = modalities{i};
        t = table(iid_keep(:), eid_keep(:), 'VariableNames', {'iid', 'eid'});
        output_file = fullfile(dir_results, sprintf('iid_eid_%s.txt', modality));
        writetable(t, output_file);
        fprintf('Wrote %s (designmat=%d post-intersect=%d dropped=%d)\n', ...
            output_file, numel(keys), height(t), numel(keys) - height(t));
    end
end

fprintf('extract_iid_eid: done (post-FEMA-intersect lists).\n');

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
