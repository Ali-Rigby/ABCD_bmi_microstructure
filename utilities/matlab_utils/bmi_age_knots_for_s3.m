function [knots, qc] = bmi_age_knots_for_s3(paths)
%BMI_AGE_KNOTS_FOR_S3 Age spline knot vector for s3_create_basis_vars_pheno.
%
% Resolution order:
%   1. BMI_AGE_KNOTS              explicit numeric vector (whitespace-separated)
%   2. BMI_AGE_KNOTS_PER_SESSION_MEDIAN=1
%                               median ab_g_dyn__visit_age per ABCD wave in the
%                               profile imaging base (ses-00A … ses-06A)
%   3. default fixed knots [9.96 12.0 14.1 16.1] (original full-sample medians)

explicit = strtrim(getenv('BMI_AGE_KNOTS'));
if ~isempty(explicit)
    knots = sscanf(explicit, '%f')';
    if isempty(knots)
        error('bmi_age_knots_for_s3:CouldNotParse', 'BMI_AGE_KNOTS: %s', explicit);
    end
    qc = struct('source', 'env_BMI_AGE_KNOTS', 'knots', knots);
    return;
end

if bmi_env_truthy('BMI_AGE_KNOTS_PER_SESSION_MEDIAN')
    [knots, qc] = bmi_age_knots_per_session_median(paths);
    return;
end

knots = [9.96 12.0 14.1 16.1];
qc = struct('source', 'default_fixed', 'knots', knots);
end

function [knots, qc] = bmi_age_knots_per_session_median(paths)
sessions = {'ses-00A', 'ses-02A', 'ses-04A', 'ses-06A'};

base_name = strtrim(getenv('BMI_IMAGING_BASE_DATAFILE_IMAGING'));
if isempty(base_name)
    base_name = 'datafile6.0_before_bfs_anthroqc_bmip_agediscrepancy_dmri_incl.txt';
end
base_file = fullfile(paths.shared_datafiles, base_name);
if ~isfile(base_file)
    error('bmi_age_knots_per_session_median:BaseFileNotFound', '%s', base_file);
end

T = readtable(base_file, 'FileType', 'text', 'Delimiter', '\t');
required = {'session_id', 'ab_g_dyn__visit_age'};
missing = required(~ismember(required, T.Properties.VariableNames));
if ~isempty(missing)
    error('bmi_age_knots_per_session_median:MissingColumns', '%s in %s', ...
        strjoin(missing, ', '), base_file);
end

age = T.ab_g_dyn__visit_age;
sid = string(T.session_id);
knots = nan(1, numel(sessions));
n_per = zeros(1, numel(sessions));

for i = 1:numel(sessions)
    mask = sid == sessions{i} & ~isnan(age);
    n_per(i) = sum(mask);
    if n_per(i) == 0
        error('bmi_age_knots_per_session_median:NoRows', 'No visit_age rows for %s', sessions{i});
    end
    knots(i) = median(age(mask));
end

qc = struct();
qc.source = 'per_session_median';
qc.base_file = base_file;
qc.sessions = sessions;
qc.n_per_session = n_per;
qc.medians = knots;

qc_dir = fullfile(paths.analysis_base, 'qc');
if ~exist(qc_dir, 'dir')
    mkdir(qc_dir);
end
qc_file = fullfile(qc_dir, 'age_knots_per_session_median.tsv');
fid = fopen(qc_file, 'w');
fprintf(fid, 'session_id\tn_visits\tmedian_visit_age\n');
for i = 1:numel(sessions)
    fprintf(fid, '%s\t%d\t%.8f\n', sessions{i}, n_per(i), knots(i));
end
fclose(fid);
qc.qc_file = qc_file;

fprintf('Age knots from per-session medians (%s):\n', base_file);
for i = 1:numel(sessions)
    fprintf('  %s  n=%d  median=%.6f\n', sessions{i}, n_per(i), knots(i));
end
fprintf('  QC: %s\n', qc_file);
end

function tf = bmi_env_truthy(varName)
tf = ismember(lower(strtrim(getenv(varName))), {'1', 'true', 'yes', 'on'});
end
