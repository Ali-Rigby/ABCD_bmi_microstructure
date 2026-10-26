function apply_s8_analysis_profile(profile, results_run)
%APPLY_S8_ANALYSIS_PROFILE Set BMI_IMAGING_* env vars for s8 / showVol runs.
%
% Profiles used in this repository:
%   bmip, perc_median
%     -> main analysis under $CEPH/<pheno>/
%   ve11_nove11_<VE11|noVE11>_<bmip|perc_median>
%     -> matched VE11 vs noVE11 sensitivity arms

if nargin >= 2 && ~isempty(results_run)
    setenv('BMI_IMAGING_RESULTS_RUN', results_run);
    setenv('BMI_IMAGING_AUTO_RUN_DATE', '0');
end

if nargin < 1 || isempty(profile)
    return;
end

setenv('BMI_IMAGING_RESULTS_ROOT', '');
setenv('BMI_IMAGING_SHARED_DATAFILES_DIR', '');
setenv('BMI_IMAGING_DATAFILES_DIR', '');
setenv('BMI_IMAGING_ABCD_SYNC', '');
setenv('BMI_IMAGING_VOXELWISE_RSI_DIR', '');

switch profile
    case 'bmip'
        setenv('BMI_PHENO_VARIABLE', 'bmip');
        setenv('BMI_IMAGING_ANALYSIS_DIR', 'bmip');
    case 'perc_median'
        setenv('BMI_PHENO_VARIABLE', 'perc_median');
        setenv('BMI_IMAGING_ANALYSIS_DIR', 'perc_median');
    otherwise
        tok_ve11_nove11 = regexp(profile, ...
            '^ve11_nove11_(VE11|noVE11)_(bmip|perc_median)$', ...
            'tokens', 'once');
        if isempty(tok_ve11_nove11)
            error(['Unknown analysis_profile=%s. Use bmip, perc_median, or ', ...
                've11_nove11_<VE11|noVE11>_<bmip|perc_median>.'], profile);
        end
        arm = tok_ve11_nove11{1};
        pheno = tok_ve11_nove11{2};
        root = 've11_vs_noVE11_6_0';
        setenv('BMI_PHENO_VARIABLE', pheno);
        setenv('BMI_SCANNER_ARM', arm);
        setenv('BMI_VE11_NOVE11_ANALYSIS_ROOT', root);
        setenv('BMI_IMAGING_ANALYSIS_DIR', sprintf('%s/%s/%s', root, pheno, arm));
        setenv('BMI_IMAGING_DATAFILES_DIR', sprintf('%s/%s/datafiles', root, pheno));
        setenv('BMI_IMAGING_SHARED_DATAFILES_DIR', sprintf('%s/shared/datafiles', root));
        setenv('BMI_IMAGING_ABCD_SYNC', ...
            '/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0');
        setenv('BMI_IMAGING_VOXELWISE_RSI_DIR', '');
        if strcmp(pheno, 'perc_median') && isempty(getenv('BMI_IMAGING_PHENO_GRID_STEP'))
            setenv('BMI_IMAGING_PHENO_GRID_STEP', '1');
        end
end
end
