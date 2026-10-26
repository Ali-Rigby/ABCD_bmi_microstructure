function R = fema_load_results(varargin)
%FEMA_LOAD_RESULTS Load FEMA estimates + mask for downstream s7/s8/s9.
%
% Name-value (one of two modes):
%   'ResultsFile', path       — estimates .mat path
%   'DirnameOut', dir         — with 'Datatype' + 'FstemImaging'
%   'Datatype', 'voxel'       — from s6 loop or filename parse
%   'FstemImaging', 'rni'     — FEMA fstem_imaging (= modality tag in filenames)
%   'DirnameImaging', path    — optional; inferred via fema_imaging_dirname if omitted
%
% Mask: loaded from results file if present (legacy); otherwise vol_info.mat
%       via fema_load_mask (same source as FEMA_process_data for voxel data).
%
% R fields: beta_hat, coeffCovar, colnames_model, mask, datatype, fstem_imaging,
%          modality, fname_results, dirname_imaging

p = inputParser;
addParameter(p, 'ResultsFile', '', @ischar);
addParameter(p, 'DirnameOut', '', @ischar);
addParameter(p, 'Datatype', '', @ischar);
addParameter(p, 'FstemImaging', '', @ischar);
addParameter(p, 'DirnameImaging', '', @ischar);
addParameter(p, 'AbcdSyncPath', '', @ischar);
parse(p, varargin{:});

if ~isempty(p.Results.ResultsFile)
    fname_results = p.Results.ResultsFile;
    [dt_parse, fs_parse] = fema_parse_results_stem(fname_results);
    datatype = pick_char(p.Results.Datatype, dt_parse);
    fstem_imaging = pick_char(p.Results.FstemImaging, fs_parse);
elseif ~isempty(p.Results.DirnameOut)
    datatype = p.Results.Datatype;
    fstem_imaging = p.Results.FstemImaging;
    if isempty(datatype) || isempty(fstem_imaging)
        error('fema_load_results: DirnameOut mode requires Datatype and FstemImaging.');
    end
    fname_results = fema_results_path(p.Results.DirnameOut, datatype, fstem_imaging);
else
    error('fema_load_results: specify ResultsFile or DirnameOut + Datatype + FstemImaging.');
end

if isempty(datatype) || isempty(fstem_imaging)
    error('fema_load_results: could not determine datatype/fstem for %s.', fname_results);
end

if ~isempty(p.Results.DirnameImaging)
    dirname_imaging = p.Results.DirnameImaging;
else
    abcd = p.Results.AbcdSyncPath;
    if isempty(abcd)
        abcd = fema_default_abcd_sync_path();
    end
    dirname_imaging = fema_imaging_dirname(datatype, fstem_imaging, abcd);
end

if ~isfile(fname_results)
    error('fema_load_results: file not found: %s', fname_results);
end

est_vars = whos('-file', fname_results);
est_names = {est_vars.name};
need = {'beta_hat', 'coeffCovar', 'colnames_model'};
missing = setdiff(need, est_names);
if ~isempty(missing)
    error('fema_load_results: %s missing: %s', fname_results, strjoin(missing, ', '));
end

S = load(fname_results, need{:});

if ismember('mask', est_names)
    Sm = load(fname_results, 'mask');
    mask = Sm.mask;
    mask_source = 'results_file';
else
    mask = fema_load_mask(dirname_imaging, datatype);
    mask_source = 'vol_info.mat';
end

R = struct();
R.beta_hat = S.beta_hat;
R.coeffCovar = S.coeffCovar;
R.colnames_model = S.colnames_model;
R.mask = mask;
R.mask_source = mask_source;
R.datatype = char(datatype);
R.fstem_imaging = char(fstem_imaging);
R.modality = R.fstem_imaging;
R.fname_results = fname_results;
R.dirname_imaging = dirname_imaging;

fprintf(['fema_load_results: %s | datatype=%s modality=%s | mask from %s | ', ...
    'imaging=%s\n'], fname_results, R.datatype, R.modality, mask_source, dirname_imaging);
end

function v = pick_char(primary, fallback)
if ~isempty(primary)
    v = char(primary);
elseif ~isempty(fallback)
    v = char(fallback);
else
    v = '';
end
end
