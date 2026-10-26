% Run FEMA for bmi analysis of 6.1 with svd of bmi basis functions

paths = get_project_paths();
fprintf('s6: analysis=%s\n', paths.analysis_base);
fprintf('s6: results_dir=%s\n', paths.results_dir);
fprintf('s6: designmats=%s\n', paths.designmats_dir);

% Specify where to store results (override with BMI_IMAGING_* env vars; see get_project_paths.m)
dirname_out = paths.results_dir;

if ~exist(dirname_out, 'dir')
    mkdir(dirname_out);
    fprintf('s6: created results directory: %s\n', dirname_out);
end

% start diary
diary_path=strcat(dirname_out,'/','diary_', datestr(now, 'yyyy-mm-dd_HHMM'));
diary(diary_path);

% Add cmig_tools to path
addpath(genpath(paths.cmig_tools))
addpath(genpath(paths.matlab_utils))

% Get path to local copy of abcd-sync using abcdConfig
% abcd_sync_path=fullfile('/space/ceph/1/ABCD/data/release_archive/6_0');
abcd_sync_path = fullfile('/space/ceph/1/ABCD/users/parekh/2025-11-19_ABCD_concat_uncompressed/6_0');

%% INPUTS TO FEMA_wrapper.m
% dirname_tabulated = fullfile(abcd_sync_path,'tabulated/img'); % directory to tabulated imaging data on abcd-sync 
% dirname_phenofiles = '/space/ceph/1/ABCD/users/arigby/2025_bmi_manuscript/phenofiles';

%To run multiple design matrices with same imaging data populate each row with path to each design matrix
designmat_dir = paths.designmats_dir;
designmat_file = dir(sprintf('%s/*.txt', designmat_dir));
designmat_file = {designmat_file.name}';
fname_design = strcat(designmat_dir, '/', designmat_file);

outdir_file = strrep(designmat_file, '.txt', '');
dirname_out = strcat(dirname_out,'/',outdir_file);

%% **TOGGLES FOR ANALYSIS (Set to 1 to enable, 0 to disable)**
do_smri_vertex = 0;  % Structural MRI Vertex-based analysis
do_dmri_vertex = 0;  % Diffusion MRI Vertex-based analysis
do_smri_voxel = 0;   % Structural MRI Voxel-based analysis
do_dmri_voxel = 1;   % Diffusion MRI Voxel-based analysis
do_external = 0;     % External analysis (ROI/behavioral data)

ranknorm = 0;

% output = 'nifti'; % toggling output format - default is 'mat'
output = 'mat';
if ranknorm
    transform_y = 'ranknorm';
else
    transform_y = 'none';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% **LOAD IMAGING DATA PATHS**

% **sMRI Vertexwise Data**
smri_vertex_modality = {'area_sm1000', 'thk_sm1000', 'sulc_sm1000', 't1_gwc_sm1000'};
smri_vertex_dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/vertexwise/smri');

% **dMRI Vertexwise Data**
dmri_vertex_modality = {'rnt_gm_sm1000', 'rnt_wm_sm1000', 'rni_gm_sm1000', 'rni_wm_sm1000', 'rnd_gm_sm1000', 'rnd_wm_sm1000', 'hnt_gm_sm1000', 'hnt_wm_sm1000', 'fni_gm_sm1000', 'fni_wm_sm1000'};
dmri_vertex_dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/vertexwise/rsi');
% **sMRI Voxelwise Data**
smri_voxel_modality = {'b0'};
smri_voxel_dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/voxelwise/smri');

% **dMRI Voxelwise Data**
dmri_voxel_modality = {'rni'};
dmri_voxel_dirname_imaging = fullfile(abcd_sync_path, 'concatenated/imaging/voxelwise/rsi');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% **VERTEXWISE ANALYSES**

% **sMRI Vertex**
if do_smri_vertex
    datatype = 'vertex';
    for m = 1:length(smri_vertex_modality)
        fstem_imaging = smri_vertex_modality{m};
        fprintf('Running FEMA for sMRI Vertex: %s\n', fstem_imaging);
        for d = 1:numel(fname_design)
            est_file = fema_run_wrapper(fstem_imaging, fname_design{d}, dirname_out{d}, ...
                smri_vertex_dirname_imaging, datatype, transform_y, output);
            fprintf('FEMA saved: %s\n', est_file);
        end
    end
end

% **dMRI Vertex**
if do_dmri_vertex
    datatype = 'vertex';
    for m = 1:length(dmri_vertex_modality)
        fstem_imaging = dmri_vertex_modality{m};
        fprintf('Running FEMA for dMRI Vertex: %s\n', fstem_imaging);
        for d = 1:numel(fname_design)
            est_file = fema_run_wrapper(fstem_imaging, fname_design{d}, dirname_out{d}, ...
                dmri_vertex_dirname_imaging, datatype, transform_y, output);
            fprintf('FEMA saved: %s\n', est_file);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% **VOXELWISE ANALYSES**

% **sMRI Voxel**
if do_smri_voxel
    datatype = 'voxel';
    for m = 1:length(smri_voxel_modality)
        fstem_imaging = smri_voxel_modality{m};
        fprintf('Running FEMA for sMRI Voxel: %s\n', fstem_imaging);
        for d = 1:numel(fname_design)
            est_file = fema_run_wrapper(fstem_imaging, fname_design{d}, dirname_out{d}, ...
                smri_voxel_dirname_imaging, datatype, transform_y, output);
            fprintf('FEMA saved: %s\n', est_file);
        end
    end
end

% **dMRI Voxel**
if do_dmri_voxel
    datatype = 'voxel';
    for m = 1:length(dmri_voxel_modality)
        fstem_imaging = dmri_voxel_modality{m};
        fprintf('Running FEMA for dMRI Voxel: %s\n', fstem_imaging);
        for d = 1:numel(fname_design)
            est_file = fema_run_wrapper(fstem_imaging, fname_design{d}, dirname_out{d}, ...
                dmri_voxel_dirname_imaging, datatype, transform_y, output);
            fprintf('FEMA saved: %s\n', est_file);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EXTERNAL ANALYSIS

if do_external

      datatype='external'; % can use txt with columns of ANY data type (e.g. ROIs, behavior) - runs mass univaraite LME across every column

      imaging_dir = dirname_phenofiles; % for using cleaned data rather than data from abcd-sync
      % modality = {'dti' 'rsi' 'smr'};
      modality = {'smr'}; % just structural measures for manuscript

      for m=1:length(modality)
            % dirname_imaging = dir(sprintf('%s/mri_y_%s_*.csv', imaging_dir, modality{m}));
            dirname_imaging = dir(sprintf('%s/*.csv', imaging_dir)); % for running all files in imaging_dir
            imaging_file = {dirname_imaging.name}';
      
            for i=1:length(imaging_file)
                  imaging_path = strcat(imaging_dir,'/',imaging_file{i});
                  fstem_imaging = strrep(imaging_file{i},'.csv','');
                  try
                        % Run FEMA
                        for d = 1:numel(fname_design)
                            est_file = fema_run_wrapper(fstem_imaging, fname_design{d}, ...
                                dirname_out{d}, imaging_path, datatype, transform_y, output);
                            fprintf('FEMA saved: %s\n', est_file);
                        end
                  catch ME
                        fprintf('Error running fema_wrapper using %s: %s', fstem_imaging, ME.message);
                        continue; % Jump to next iteration
                  end
            end 
      end
end