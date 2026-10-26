function paths = get_project_paths()
%GET_PROJECT_PATHS Unified BMI imaging pipeline (bmip / perc_median).
%
% Shared tables: paths.shared_datafiles  (<ceph>/shared/datafiles)
% Profile I/O:   paths.datafiles         (<ceph>/<pheno>/datafiles)


this_file_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_file_dir);

pheno = getenv('BMI_PHENO_VARIABLE');
if isempty(pheno)
    error('Set BMI_PHENO_VARIABLE to bmip or perc_median.');
end

if isempty(getenv('BMI_IMAGING_ANALYSIS_DIR'))
    switch pheno
        case 'bmip'
            setenv('BMI_IMAGING_ANALYSIS_DIR', 'bmip');
        case 'perc_median'
            setenv('BMI_IMAGING_ANALYSIS_DIR', 'perc_median');
        otherwise
            error('BMI_PHENO_VARIABLE must be bmip or perc_median.');
    end
end

ensure_run_tag_env();

cephRoot = getenv('BMI_IMAGING_CEPH_ROOT');
if isempty(cephRoot)
    cephRoot = '/space/ceph/1/ABCD/users/arigby/2026_bmi_imaging_analysis';
end

analysisDir = getenv('BMI_IMAGING_ANALYSIS_DIR');
if analysisDir(1) == '/'
    analysisBase = analysisDir;
else
    analysisBase = fullfile(cephRoot, analysisDir);
end

paths = struct();
paths.project_root = repo_root;
paths.matlab_utils = fullfile(repo_root, 'utilities', 'matlab_utils');
paths.cmig_tools = '/usr/pubsw/packages/cmig/cmig_tools';
paths.ceph_project = cephRoot;
sharedDir = getenv('BMI_IMAGING_SHARED_DATAFILES_DIR');
if isempty(sharedDir)
    paths.shared_datafiles = fullfile(cephRoot, 'shared', 'datafiles');
elseif sharedDir(1) == '/'
    paths.shared_datafiles = sharedDir;
else
    paths.shared_datafiles = fullfile(cephRoot, sharedDir);
end
paths.analysis_base = analysisBase;
paths.pheno_variable = pheno;

% showVol atlas underlay (override with BMI_IMAGING_SHOWVOL_ATLAS)
showvolAtlas = getenv('BMI_IMAGING_SHOWVOL_ATLAS');
if isempty(showvolAtlas)
    showvolAtlas = '/home/dpecheva/tools/internal/showVolData/Atlas/6.0_ABCD3_cor10';
end
paths.showvol_atlas = showvolAtlas;

resultsRoot = getenv('BMI_IMAGING_RESULTS_ROOT');
if isempty(resultsRoot)
    paths.results_root = fullfile(analysisBase, 'results');
elseif resultsRoot(1) == '/'
    paths.results_root = resultsRoot;
else
    paths.results_root = fullfile(cephRoot, resultsRoot);
end

resultsRun = getenv('BMI_IMAGING_RESULTS_RUN');
if isempty(resultsRun)
    resultsRun = ['results_' run_tag()];
    setenv('BMI_IMAGING_RESULTS_RUN', resultsRun);
end
paths.results_run = resultsRun;
paths.results_dir = fullfile(paths.results_root, resultsRun);

datafilesDir = getenv('BMI_IMAGING_DATAFILES_DIR');
if isempty(datafilesDir)
    paths.datafiles = fullfile(analysisBase, 'datafiles');
elseif datafilesDir(1) == '/'
    paths.datafiles = datafilesDir;
else
    paths.datafiles = fullfile(cephRoot, datafilesDir);
end
paths.profile_datafiles = paths.datafiles;

designmatsDir = getenv('BMI_IMAGING_DESIGNMATS_DIR');
designmatSub = getenv('BMI_IMAGING_DESIGNMATS_SUBDIR');
if isempty(designmatsDir)
    if isempty(designmatSub)
        paths.designmats_dir = fullfile(analysisBase, 'designmats');
    else
        paths.designmats_dir = fullfile(analysisBase, 'designmats', designmatSub);
    end
elseif designmatsDir(1) == '/'
    paths.designmats_dir = designmatsDir;
else
    paths.designmats_dir = fullfile(cephRoot, designmatsDir);
end

switch pheno
    case 'bmip'
        paths.pheno_linear_term = 'bmip';
        paths.pheno_xname = 'bmi';
        paths.pheno_vec_file = 'bmi_percentile_vector.txt';
        paths.pheno_bf_table = 'table_bf_bmi_perc.txt';
        paths.pheno_bf_grep = 'bf_bmi';
    case 'perc_median'
        paths.pheno_linear_term = 'perc_median';
        paths.pheno_xname = 'perc_median';
        paths.pheno_vec_file = 'perc_median_vector.txt';
        paths.pheno_bf_table = 'table_bf_perc_median.txt';
        paths.pheno_bf_grep = 'bf_perc_median';
    otherwise
        error('BMI_PHENO_VARIABLE must be bmip or perc_median.');
end

paths.run_tag = getenv('BMI_IMAGING_RUN_TAG');

end

function ensure_run_tag_env()
auto = getenv('BMI_IMAGING_AUTO_RUN_DATE');
if isempty(auto)
    auto = '1';
end
if ~ismember(lower(auto), {'1', 'true', 'yes', 'on'})
    return;
end
if ~isempty(getenv('BMI_IMAGING_DESIGNMATS_DIR'))
    return;
end
tag = run_tag();
if isempty(getenv('BMI_IMAGING_DESIGNMATS_SUBDIR'))
    setenv('BMI_IMAGING_DESIGNMATS_SUBDIR', ['designmats_' tag]);
end
if isempty(getenv('BMI_IMAGING_RESULTS_RUN'))
    setenv('BMI_IMAGING_RESULTS_RUN', ['results_' tag]);
end
end

function tag = run_tag()
tag = strtrim(getenv('BMI_IMAGING_RUN_TAG'));
if isempty(tag)
    tag = datestr(datetime('today'), 'yyyy_mm_dd');
    setenv('BMI_IMAGING_RUN_TAG', tag);
end
end
