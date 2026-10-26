function runs = fema_list_results(dirname_out, varargin)
%FEMA_LIST_RESULTS Discover FEMA estimate files under a designmat results folder.
%
% Optional name-value (for bare FEMA_estimates.mat with no parseable prefix):
%   'DefaultDatatype'   — e.g. 'voxel'
%   'DefaultFstem'      — e.g. 'rni' (same as fstem_imaging / modality tag)

p = inputParser;
addParameter(p, 'DefaultDatatype', '', @ischar);
addParameter(p, 'DefaultFstem', '', @ischar);
parse(p, varargin{:});
default_datatype = p.Results.DefaultDatatype;
default_fstem = p.Results.DefaultFstem;

runs = struct('fname_results', {}, 'datatype', {}, 'fstem_imaging', {}, 'modality', {});
seen = {};

add_run = @(fname, dt, fs) assign_run(fname, dt, fs);

% Prefixed FEMA 4.x estimates
est_files = dir(fullfile(dirname_out, 'FEMA_wrapper_output_*_estimates.mat'));
for i = 1:numel(est_files)
    add_run(est_files(i).name, [], []);
end

% Legacy single-file outputs (exclude Wald / design sidecars)
leg_files = dir(fullfile(dirname_out, 'FEMA_wrapper_output_*.mat'));
for i = 1:numel(leg_files)
    name = leg_files(i).name;
    if contains(name, '_estimates.mat') || contains(name, '_Wald_') || contains(name, 'designMatrix')
        continue;
    end
    [dt_chk, fs_chk] = fema_parse_results_stem(name);
    if ~isempty(dt_chk) && ~isempty(fs_chk) && isfile(fullfile(dirname_out, ...
            sprintf('FEMA_wrapper_output_%s_%s_estimates.mat', dt_chk, fs_chk)))
        continue;
    end
    add_run(name, [], []);
end

% Bare FEMA_estimates.mat (current 4.x without outPrefix)
if isfile(fullfile(dirname_out, 'FEMA_estimates.mat'))
    add_run('FEMA_estimates.mat', default_datatype, default_fstem);
end

    function assign_run(fname, dt_in, fs_in)
        if ismember(fname, seen)
            return;
        end
        seen{end + 1} = fname; %#ok<AGROW>
        if isempty(dt_in)
            [dt, fs] = fema_parse_results_stem(fname);
        else
            dt = dt_in;
            fs = fs_in;
        end
        if strcmp(fname, 'FEMA_estimates.mat') && (isempty(dt) || isempty(fs))
            if isempty(default_datatype) || isempty(default_fstem)
                warning(['fema_list_results: skipping %s (set DefaultDatatype/DefaultFstem).'], ...
                    fullfile(dirname_out, fname));
                return;
            end
            dt = default_datatype;
            fs = default_fstem;
        end
        if isempty(dt)
            return;
        end
        r = struct();
        r.fname_results = fullfile(dirname_out, fname);
        r.datatype = dt;
        r.fstem_imaging = fs;
        r.modality = fs;
        runs(end + 1) = r; %#ok<AGROW>
    end
end
