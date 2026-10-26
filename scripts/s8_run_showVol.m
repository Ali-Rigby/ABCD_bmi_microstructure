% showVol for unified BMI imaging pipeline (perc_median / bmip splines).

pheno_profile = 'perc_median';              % 'bmip' | 'perc_median'
results_run   = 'results_2026_06_29'; % e.g. 'results_2026_06_29'; '' = env / auto date

if ~isempty(pheno_profile)
    setenv('BMI_PHENO_VARIABLE', pheno_profile);
    setenv('BMI_IMAGING_ANALYSIS_DIR', pheno_profile);
end
if ~isempty(results_run)
    setenv('BMI_IMAGING_RESULTS_RUN', results_run);
end

paths = get_project_paths();
pheno_xname = paths.pheno_xname;
pheno_linear = paths.pheno_linear_term;
fprintf('s8 (pheno): %s\n', pheno_xname);
fprintf('Results dir: %s\n', paths.results_dir);

%% Design matrix and modality settings
desmatnum = '1f'; % which designmat_file you want, e.g. '8m'
modality = 'rni'; % FEMA fstem_imaging (= modality tag in output filenames)
datatype = 'voxel'; % voxel | vertex — must match s6

%% path to FEMA results directory (override with BMI_IMAGING_RESULTS_RUN or results_run above)
dirname_out = paths.results_dir;
if ~exist(dirname_out, 'dir')
    error('Results directory not found: %s\nSet results_run at top of script or BMI_IMAGING_RESULTS_RUN in env.', dirname_out);
end

designmat_file = list_s8_designmat_dirs(dirname_out);
if isempty(designmat_file)
    error('No designmat* folders under %s', dirname_out);
end
fprintf('Designmat folders found: %s\n', strjoin(designmat_file, ', '));
dirname_results = fullfile(dirname_out, designmat_file);

%% Add cmig_tools and utility_functions to path
addpath(genpath(paths.cmig_tools));
addpath(genpath(paths.matlab_utils));

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

%% Define specific pwaldlims for each xname for Wald stats (tune pheno range to your BF grid)
% rni — pheno = % distance from CDC median
pwaldlims_dict = struct(pheno_xname, [0, 30], 'age', [0, 150], 'pds', [0, 5]);
pwaldlims_Bon_dict = pwaldlims_dict;
pwaldlims_FDR_dict = pwaldlims_dict;

% t2 (uncomment and tune if needed)
% pwaldlims_dict = struct(pheno_xname, [-50, 100], 'age', [0, 50], 'pds', [0, 5]);
% pwaldlims_Bon_dict = pwaldlims_dict;
% pwaldlims_FDR_dict = pwaldlims_dict;

%pwaldlims_dict = struct('bmi', [0, 30]); % raw rni
%pwaldlims_Bon_dict = struct('bmi', [0, 30]); % raw rni
%pwaldlims_FDR_dict = struct('bmi', [0, 30]); % raw rni

%% Voxel limits (toggle on / off) % rni
limits = [-1 1]; % colormap limits
phenolim = [-0.0008 0.0008];
agelim = [-0.01 0.01]; 
pdslim = [-0.008 0.008]; 
sexlim = [-0.01, 0.01];

%% Voxel limits (toggle on / off) % t2
%limits = [-1 1]; % colormap limits
%phenolim = [-0.8 0.8]; 
%agelim = [-10.0 10.0]; 
%pdslim = [-5.0 5.0]; 
%sexlim = [-0.01, 0.01];

%% FEMA settings
showPval = false;               % If showPval==true, creates a map of log10 p-values
cmap_linear = redblackblue;     % Colormap
cmap_splines = redblackblue;    % Colormap for splines
cmap_wald = cmap_linear(200:end, :); % Use a subset of the colormap
bg_vol = [];                    % Set background image to overlay statistics. Default: atlas T1
CSFprob = [];                   % Probabilistic threshold for masking CSF e.g. if CSFprob=0.8, only voxels in which 80% of participants labelled that voxel as CSF will be masked (Default: no masking)
interpMethod = [];              % By default uses linear interpolation to go from 2mm to 1mm voxels. Specify 'nearest' to show actual resolution
atlasVersion = 'ABCD3';

%% Load basis function tables

% Determine if the designmat is for males or females based on desmatnum
dirname_basis = paths.datafiles;
shared_basis = paths.shared_datafiles;

tbl_bf_pheno = readtable(fullfile(dirname_basis, paths.pheno_bf_table));
tbl_bf_age = readtable(fullfile(shared_basis, 'table_bf_age.txt'));
tbl_bf_pds = readtable(fullfile(shared_basis, 'table_bf_pds.txt'));

%% Load vectors of variables of interest (match s3 BF grid extent)
vec_file = fullfile(dirname_basis, paths.pheno_vec_file);
val_table = readtable(vec_file, 'Delimiter', '\t');
phenovec = pheno_bf_xvec(tbl_bf_pheno, paths);
if height(tbl_bf_pheno) ~= numel(phenovec)
    error('BF table rows (%d) != phenovec length (%d). Re-run s3.', height(tbl_bf_pheno), numel(phenovec));
end
fprintf('Pheno eval grid: n=%d, [%g : %g : %g]\n', numel(phenovec), phenovec(1), phenovec(2)-phenovec(1), phenovec(end));
agevec = (8.3:0.1:17.8);
pdsvec = linspace(1, 5, 101);

% Spline slice indices (must match BF grid length for active pheno)
switch paths.pheno_variable
    case 'bmip'
        pheno_slice_list = [2, 11:10:91, 100];
    case 'perc_median'
        pheno_slice_list = [1, 3:5:min(163, numel(phenovec)), numel(phenovec)];
    otherwise
        pheno_slice_list = [2, round(linspace(2, numel(phenovec)-1, 12)), numel(phenovec)];
end
pheno_slice_list = pheno_slice_list(pheno_slice_list >= 1 & pheno_slice_list <= numel(phenovec));
pheno_slice_list = unique(pheno_slice_list);

%% loop through designmats
abcd_sync_path = fema_default_abcd_sync_path();
dirname_imaging = fema_imaging_dirname(datatype, modality, abcd_sync_path);

des_idx = find(contains(designmat_file, sprintf('designmat%s_', desmatnum)))';
if isempty(des_idx)
    fprintf('Available designmat folders:\n');
    fprintf('  %s\n', designmat_file{:});
    error('No folder matches desmatnum=%s under %s', desmatnum, dirname_out);
end

for des = des_idx

    est_file = fema_results_path(dirname_results{des}, datatype, modality);
    R = fema_load_results('ResultsFile', est_file, ...
        'Datatype', datatype, 'FstemImaging', modality, 'DirnameImaging', dirname_imaging);

    beta_hat = R.beta_hat;
    coeffCovar = R.coeffCovar;
    mask = R.mask;
    colnames_model = R.colnames_model;
    fname_results = R.fname_results;

    fprintf('s8: loaded %s (datatype=%s modality=%s)\n', est_file, R.datatype, R.modality);

        % Indices for variables of interest
        pheno_idx = find(strcmp(colnames_model, pheno_linear));
        age_idx = find(strcmp(colnames_model, 'ab_g_dyn__design_age__event'));
        pds_idx = find(contains(colnames_model, 'pds_categ_y_p_all'));
        sex_idx = find(contains(colnames_model, 'sexF'));
        Spheno_idx = find(find_cell(regexp(colnames_model, ['^', paths.pheno_bf_grep], 'match')));
        Sage_idx = find(find_cell(regexp(colnames_model,'^bf_age','match')));
        Spds_idx = find(find_cell(regexp(colnames_model,'^bf_pds','match')));
        %Sbmixage_idx = find(find_cell(regexp(colnames_model,'bf_bmi.*_interview_age','match')));
        %Sbmixpds_idx = find(find_cell(regexp(colnames_model,'bf_bmi.*pds','match')));
        intercept_idx = find(strcmp(colnames_model, 'intercept'));
        
        % For Wald stats
        %xnames = {pheno_xname, 'age', 'pds'};
        xnames = {pheno_xname};

        if strcmp(datatype,'voxel')
            clear vols;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% WALD STAT MAPS

            % Loop through Wald stats
            vols = []; % Initialize vols for voxelwise stats
            for i = 1:length(xnames)
                xname = xnames{i};
                pwaldlimits = pwaldlims_dict.(xname);
                pwaldlimits_Bon = pwaldlims_Bon_dict.(xname);
                pwaldlimits_FDR = pwaldlims_FDR_dict.(xname);

                % Process Wald stats
                vols = process_wald_stats(vols, dirname_results{des}, modality, xname, designmat_file{des}, pwaldlimits, pwaldlimits_Bon, pwaldlimits_FDR, showPval, cmap_wald, bg_vol, CSFprob, interpMethod, atlasVersion);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % LINEAR MAPS

            % idx_list = {pheno_idx, age_idx, sex_idx, pds_idx};
            % limits_list = {phenolim, agelim, sexlim, pdslim};
            % cmap = redblackblue;

            % vols = process_linear_maps(vols, vol_beta_hat, modality, designmat_file{des}, colnames_model, idx_list, limits_list, showPval, cmap_linear, bg_vol, CSFprob, interpMethod, atlasVersion);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SPLINE MAPS
            
            if any(contains(colnames_model, 'bf_'))
                %spline_vars = {pheno_xname, 'age', 'pds'};
                spline_vars = {pheno_xname};
                tbl_list = {tbl_bf_pheno, tbl_bf_age, tbl_bf_pds};
                vec_list = {phenovec, agevec, pdsvec};
                idx_list = {Spheno_idx, Sage_idx, Spds_idx};
                limits_list = {phenolim, agelim, pdslim};
                slice_list = {pheno_slice_list, [3, 8:10:88, 96], [1, 13, 26, 38, 51, 63, 76, 88, 101]};
                %slice_list = {[2, 11:10:81, 86, 91, 96, 100]};
                
                for i = 1:length(spline_vars)
                    if isempty(idx_list{i})
                        continue;
                    end
                    vols = process_spline_slices(vols, tbl_list{i}, vec_list{i}, beta_hat, idx_list{i}, ...
                                                intercept_idx, mask, coeffCovar, slice_list{i}, ...
                                                modality, designmat_file{des}, colnames_model, limits_list{i}, ...
                                                showPval, cmap_splines, bg_vol, CSFprob, interpMethod, atlasVersion, spline_vars{i});
                end
            end

            % if any(contains(colnames_model, 'bf_'))
            %     % BMI
            %     vols = process_spline_slices(vols, tbl_bf_pheno, phenovec, beta_hat, Spheno_idx, intercept_idx, mask, coeffCovar, [2, 11:10:91, 100], modality, designmat_file{des}, colnames_model, phenolim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, 'bmi');

            %     % Age
            %     vols = process_spline_slices(vols, tbl_bf_age, agevec, beta_hat, Sage_idx, intercept_idx, mask, coeffCovar, [3, 8:10:88, 97], modality, designmat_file{des}, colnames_model, agelim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, 'age');

            %     % PDS:
            %     % vols = process_spline_slices(vols, tbl_bf_pds, pdsvec, beta_hat, Spds_idx, intercept_idx, mask, coeffCovar, [1, 13, 26, 38, 51, 63, 76, 88, 101], modality, designmat_file{des}, colnames_model, pdslim, showPval, cmap, bg_vol, CSFprob, interpMethod, atlasVersion, 'pds');
            % end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% RUN SHOWVOL
            coords = [0 -17 2];

            if isempty(vols)
                warning('s8: no volumes to display for %s (check s7 Wald outputs and spline processing).', designmat_file{des});
            else
                fprintf('s8: opening showVol with %d volumes for %s\n', numel(vols), designmat_file{des});
                showVol(vols, eval(vname_cFA), eval(vname_t1), coords);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % VERTEX MAPS

        if strcmp(datatype,'vertex')

            legendPosition = 'south';
            ico = 5; % ico number
            polarity_Wald = 1;  % Wald stats require polarity = 1
            polarity_Other = 2; % Everything else requires polarity = 2

            % % Values for CT, SA, GWC
            % % Vertex limits (toggle on / off) 
            % limits = [-1 1]; % colormap limits
            % phenolim = [-.06 .06];
            % agelim = [-0.25 0.25];
            % pdslim = [-0.1 0.1];
            % sexlim = [-0.01, 0.01];

            % % Define specific pwaldlims for each xname
            % pwaldlims_dict = struct('bmi', [0, 30], 'age', [0, 300], 'pds', [0, 10]);

            % Values for RNT, RNI, RND
            % Vertex limits (toggle on / off)
            limits = [-1 1]; % colormap limits
            %phenolim = [-.015 .015]; % rn gray
            phenolim = [-.01 .01]; % rn white
            %phenolim = [-.0004 .0004] % raw
            agelim = [-0.25 0.25];
            pdslim = [-0.1 0.1];
            sexlim = [-0.01, 0.01];

            % Define specific pwaldlims for each xname
            pwaldlims_dict = struct(pheno_xname, [-100, 200], 'age', [0, 300], 'pds', [0, 10]); % gm
            %pwaldlims_dict = struct('bmi', [0, 25], 'age', [0, 300], 'pds', [0, 10]); % wm

            % % Values for HNT
            % % Vertex limits (toggle on / off)
            % limits = [-1 1]; % colormap limits
            % phenolim = [-.05 .05];
            % agelim = [-0.25 0.25];
            % pdslim = [-0.1 0.1];
            % sexlim = [-0.01, 0.01];

            % % Define specific pwaldlims for each xname
            % pwaldlims_dict = struct('bmi', [0, 15], 'age', [0, 300], 'pds', [0, 10]);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % WALD STAT MAPS

            % Loop through Wald stats
            for i = 1:length(xnames)
                xname = xnames{i};
                pwaldlimits = pwaldlims_dict.(xname);

                % Load corresponding Wald stat files
                splines_file = sprintf('%s/%s_Wald_splines_CI_output_%s.mat', dirname_results{des}, modality, xname);
                if exist(splines_file, 'file')
                    load(splines_file);
                else
                    warning('Splines output file not found for %s', xname);
                    continue;
                end

                fmax = pwaldlimits(2);
                fmid = fmax/2;
                fmin = 0; % middle value for plotting purposes
                fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function

                clim = pwaldlimits;

                %for FEMA_run_showSurf
                str = sprintf('\\nabla%s ', modality);
                cmap = colormap_blackred();
                ico = 5; % ico number
                icnum = ico + 1; % index for ico number (icnum = ico + 1)
                climMax = pwaldlimits(2);
                climMin = pwaldlimits(1);
                subplotDims = [3 2];
                figPosition = [10 10 16 15];
                viewname = {'bottom', 'top', 'left', 'right', 'right', 'left'};
                includemat = [1 1; 1 1; 1 0; 0 1; 1 0; 0 1];


                % Linear Wald Stat
                if exist('surf_logp_l', 'var')
                    vertvals = surf_logp_l;
                    statname = sprintf('%s -log10p Wald stat linear, %s', xname, designmat_file{des});
                    FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                        'colormap', colormap_blackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Wald);
                end

                % Nonlinear Wald Stat
                if exist('surf_logp_s', 'var')
                    vertvals = surf_logp_s;
                    statname = sprintf('%s -log10p Wald stat nonlinear, %s', xname, designmat_file{des});
                    % FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    %     'colormap', colormap_blackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Wald);
                    showSurfPlot(vertvals, str, cmap, polarity_Wald, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end

                % Combined Wald Stat
                if exist('surf_logp_ls', 'var')
                    vertvals = surf_logp_ls;
                    statname = sprintf('%s -log10p Wald stat linear + nonlinear, %s', xname, designmat_file{des});
                %     FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                %         'colormap', colormap_blackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Wald);
                    showSurfPlot(vertvals, str, cmap, polarity_Wald, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end

                % Bonferroni-corrected Linear Wald Stat
                if exist('surf_logp_l_Bon', 'var')
                    vertvals = surf_logp_l_Bon;
                    statname = sprintf('%s Bonferroni-corrected -log10p Wald stat linear, %s', xname, designmat_file{des});
                    showSurfPlot(vertvals, str, cmap, polarity, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end

                % Bonferroni-corrected Nonlinear Wald Stat
                if exist('surf_logp_s_Bon', 'var')
                    vertvals = surf_logp_s_Bon;
                    statname = sprintf('%s Bonferroni-corrected -log10p Wald stat nonlinear, %s', xname, designmat_file{des});
                %     FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                %         'colormap', colormap_blackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Wald);
                    showSurfPlot(vertvals, str, cmap, polarity_Wald, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end

                % Bonferroni-corrected Combined Wald Stat
                if exist('surf_logp_ls_Bon', 'var')
                    vertvals = surf_logp_ls_Bon;
                    statname = sprintf('%s Bonferroni-corrected -log10p Wald stat linear + nonlinear, %s', xname, designmat_file{des});
                    % FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    %     'colormap', colormap_blackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Wald);
                    showSurfPlot(vertvals, str, cmap, polarity_Wald, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LINEAR MAPS

            % Compute beta hat for BMI
            fmax= phenolim(2);
            fmid = fmax/2;
            fmin = 0; % middle value for plotting purposes
            fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            clim = phenolim;
            if ~isempty(pheno_idx)
                vertvals = beta_hat(pheno_idx, :);
                statname = sprintf('%s beta_hat, %s', pheno_xname, designmat_file{des});
                FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            end

            % Compute beta hat for Age
            fmax= agelim(2);
            fmid = fmax/2;
            fmin = 0; % middle value for plotting purposes
            fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            clim = agelim;
            if ~isempty(age_idx)
                vertvals = beta_hat(age_idx, :);
                statname = sprintf('age beta_hat, %s', designmat_file{des});
                FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            end

            % Compute beta hat for PDS
            fmax= pdslim(2);
            fmid = fmax/2;
            fmin = 0; % middle value for plotting purposes
            fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            clim = pdslim;
            if ~isempty(pds_idx)
                vertvals = beta_hat(pds_idx, :);
                statname = sprintf('pds beta_hat, %s', designmat_file{des});
                FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % % % SPLINE MAPS

            % % Compute Spline Values for BMI
            % fmax= phenolim(2);
            % fmid = fmax/2;
            % fmin = 0; % middle value for plotting purposes
            % fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            % clim = phenolim;

            %for FEMA_run_showSurf
            str = sprintf('\\nabla%s ', modality);
            cmap = blueblackred();
            ico = 5; % ico number
            icnum = ico + 1; % index for ico number (icnum = ico + 1)
            polarity = 2;
            climMax = phenolim(2);
            climMin = phenolim(1);
            subplotDims = [3 2];
            figPosition = [10 10 16 15];
            viewname = {'bottom', 'top', 'left', 'right', 'right', 'left'};
            includemat = [1 1; 1 1; 1 0; 0 1; 1 0; 0 1];

            if ~isempty(Spheno_idx)
                [surf_dbf, ~, ~, surf_bf, ~] = FEMA_convert_splinesurfs(tbl_bf_pheno, phenovec, beta_hat, Spheno_idx, intercept_idx, mask, coeffCovar);
                slice_indices = [2, 11:10:91, 100];
                for slice = slice_indices
                    vertvals = surf_dbf(slice, :);
                    statname = sprintf('%s ~ %s %.1f, %s', modality, pheno_xname, phenovec(slice), designmat_file{des});
                    % FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
                    %     'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
                    showSurfPlot(vertvals, str, cmap, polarity, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
                end
            end
            

            % % Compute Spline Values for BMI
            % fmax= phenolim(2);
            % fmid = fmax/2;
            % fmin = 0; % middle value for plotting purposes
            % fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            % clim = phenolim;
            % if ~isempty(Spheno_idx)
            %     [surf_dbf, ~, ~, surf_bf, ~] = FEMA_convert_splinesurfs(tbl_bf_pheno, phenovec, beta_hat, Spheno_idx, intercept_idx, mask, coeffCovar);
            %     for slice = 1:10:min(size(surf_dbf, 1), length(phenovec))
            %         vertvals = surf_dbf(slice, :);
            %         statname = sprintf('%s ~ bmi %.1f, %s', modality, phenovec(slice), designmat_file{des});
            %         FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
            %             'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            %     end
            % end

            % % Compute Spline Values for Age
            % fmax= agelim(2);
            % fmid = fmax/2;
            % fmin = 0; % middle value for plotting purposes
            % fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            % clim = agelim;

            % %for FEMA_run_showSurf
            % str = sprintf('\\nabla%s ', modality);
            % cmap = blueblackred();
            % ico = 5; % ico number
            % icnum = ico + 1; % index for ico number (icnum = ico + 1)
            % polarity = 2;
            % climMax = agelim(2);
            % climMin = agelim(1);
            % subplotDims = [3 2];
            % figPosition = [10 10 16 15];
            % viewname = {'bottom', 'top', 'left', 'right', 'right', 'left'};
            % includemat = [1 1; 1 1; 1 0; 0 1; 1 0; 0 1];

            % if ~isempty(Sage_idx)
            %     [surf_dbf, ~, ~, surf_bf, ~] = FEMA_convert_splinesurfs(tbl_bf_age, agevec, beta_hat, Sage_idx, intercept_idx, mask, coeffCovar);
            %     slice_indices = [3, 8:10:88, 97];
            %     for slice = slice_indices
            %     %for slice = 1:10:min(size(surf_dbf, 1), length(agevec))
            %         vertvals = surf_dbf(slice, :);
            %         % statname = sprintf('%s ~ age %.1f, %s', modality, agevec(slice), designmat_file{des});
            %         str = sprintf('%s ~ age %.1f, %s', modality, agevec(slice), designmat_file{des});
            %         %FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
            %          %   'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            %         showSurfPlot(vertvals, str, cmap, polarity, climMin, climMax, ico, 'viewname', viewname, 'includemat', includemat, 'subplotDims', subplotDims, 'figPosition', figPosition);
            %     end
            % end
            
            % % Compute spline values for PDS
            % fmax= pdslim(2);
            % fmid = fmax/2;
            % fmin = 0; % middle value for plotting purposes
            % fvals = [fmin fmid fmax]; % this will be passed to the SurfView_show_new function
            % clim = pdslim;
            % if ~isempty(Spds_idx)
            %     [surf_dbf, ~, ~, surf_bf, ~] = FEMA_convert_splinesurfs(tbl_bf_pds, pdsvec, beta_hat, Spds_idx, intercept_idx, mask, coeffCovar);
            %     for slice = 1:25:min(size(surf_dbf, 1), length(pdsvec))
            %         vertvals = surf_dbf(slice, :);
            %         statname = sprintf('%s ~ pds %.1f, %s', modality, pdsvec(slice), designmat_file{des});
            %         FEMA_run_showSurf(vertvals, statname, fvals, clim, ...
            %             'colormap', blueblackred(), 'legendPosition', legendPosition, 'ico', ico, 'polarity', polarity_Other);
            %     end
            % end
        end
end

function names = list_s8_designmat_dirs(dirname_out)
    d = dir(fullfile(dirname_out, 'designmat*'));
    d = d([d.isdir]);
    names = {d.name}';
    keep = false(size(names));
    for i = 1:numel(names)
        keep(i) = ~isempty(regexp(names{i}, '^designmat\d+[fm]_', 'once'));
    end
    names = sort(names(keep));
end