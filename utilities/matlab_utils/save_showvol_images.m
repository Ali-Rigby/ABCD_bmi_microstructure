function save_showvol_images(coords, outdir, fstem, varargin)
% save_showvol_images
% -------------------------------------------------------------------------
% Capture brain volume screenshots from an active showVol GUI instance with
% optional multi-parcellation ROI overlays, UI element suppression, and a
% structured directory layout.
%
% CORE BEHAVIOR:
%   • Iterates over all currently visible volumes in the showVol window.
%   • For each user-supplied coordinate (RCS or LPH), updates the slice view
%     (coronal / sagittal / axial based on current orientation) then saves a PNG.
%   • Optionally overlays selected ROIs (outline only) from one or more
%     parcellations, using shortcut keywords for bulk expansion.
%   • Suppresses GUI elements (colorbar, crosshairs, buttons, labels) to
%     produce publication-ready images if desired.
%   • Creates a nested output directory structure: 
%       <outdir>/screenshots/<DesignmatName>/<Modality>/<orientation>/<coord_string>/
%     where <coord_string> encodes the voxel/LPH coordinate.
%
% USAGE (minimal):
%   save_showvol_images([97 99 148], outdir, 'RNI')
%
% USAGE (with overlays & custom options):
%   save_showvol_images([97 99 148; 95 102 155], outdir, labelsCell, ...
%       'CoordSpace','RCS', 'drawROI', true, 'parcellation', {'aseg','fiber'}, ...
%       'roiNames', struct('aseg', {{'Left-Caudate','Right-Caudate'}}, 'fiber', {{'all_fiber'}}), ...
%       'lw', 2, 'fixedColor', [0 0.45 0.75], 'DesignmatName','designmat1f', 'Modality','rni');
%
% REQUIRED ARGUMENTS:
%   coords        (N x 3 numeric)  Coordinates in chosen space (LPH or RCS).
%   outdir        (char/string)    Root directory for screenshot output.
%   fstem         (char/string OR cellstr) Filename stem(s). If single, replicated
%                                      across all volumes. If a cell array, must
%                                      match number of visible volumes.
%
% KEY OPTIONAL NAME–VALUE PAIRS:
%   'CoordSpace'      : 'RCS' (default) or 'LPH'. Determines coordinate interpretation
%                       and filename formatting. RCS uses 3-digit zero padding; LPH uses
%                       2-digit padding with 'N' prefix for negatives.
%   'drawROI'         : logical (default false). If true, overlays selected ROIs.
%   'parcellation'    : cell array of parcellation names present in handles.anat (e.g. {'aseg','fiber'}).
%   'roiNames'        : struct mapping parcellation -> cell array of ROI names or shortcut keywords.
%   'lw'              : outline line width (default 1).
%   'fixedColor'      : [R G B] color for all outlines (default white).
%   'hideUI'          : master switch for hiding auxiliary UI components (default true).
%   'hideCrosshairs'  : hide crosshair overlays (default true).
%   'hideColorbar'    : hide colorbar (default true); set false to generate a custom left-side bar.
%   'hideROIText'     : hide ROI name text objects (default true).
%   'hideButtons'     : remove interactive buttons (default true).
%   'hideVolName'     : hide volume name label (default true).
%   'DesignmatName'   : string tag inserted into directory hierarchy (default 'designmat_unknown').
%   'Modality'        : string tag for modality directory (default 'modality_unknown').
%
% SHORTCUT ROI KEYWORDS (case-insensitive):
%   'subcortical_all', 'all_fiber', 'all_aparcaseg', 'all_thalamus', 'all_pauli'
%   These expand to all ROIs for the given parcellation (excluding index 1 if that
%   is a header / background entry).
%
% OUTPUT ARTIFACTS:
%   • PNG files per volume × coordinate: 
%       <modality>_<coord>_<fstem>_vol###.png
%     Example (RCS): rni_X097_Y099_Z148_RNI_vol001.png
%   • Nested directories grouping by design matrix, modality, orientation, and coordinate.
%
% COLORBAR NOTES:
%   If 'hideColorbar' is false and the volume has fields 'colormap' and 'limits', a custom
%   vertical colorbar is rendered at a fixed position (left side), showing only min/max ticks.
%   Existing colorbars tagged 'MyColorbar' are removed before recreation.
%
% PERFORMANCE & LIMITATIONS:
%   • Assumes an existing showVol figure tagged 'showVol'; does not instantiate new data.
%   • Outline alpha is forced (overlayAlpha = 0) to emphasize outlines without fill.
%   • ROI overlays accumulate across volumes unless the underlying handles are reset.
%   • If multiple showVol windows are open, only the first encountered is used.
%
% ERROR HANDLING / WARNINGS:
%   • Throws if no showVol window found.
%   • Warns when parcellation requested but missing, or ROI names absent.
%   • Validates fstem length versus number of visible volumes.
%
% COMMON TROUBLESHOOTING:
%   • Empty overlays: confirm parcellation present in handles.anat.* and ROI names match uiNames.
%   • Misaligned coordinates: ensure 'CoordSpace' matches the format of 'coords'.
%   • Missing colorbar: set 'hideColorbar', false and verify volume struct has 'colormap' & 'limits'.
%
% SEE ALSO: showVol, anatomyDrawRoiOutline_save_images
% -------------------------------------------------------------------------

%% --------------------
% Parse inputs
%% --------------------
p = inputParser;

% Core functionality
p.addParameter('drawROI', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('parcellation', {}, @(x) iscell(x) || isstring(x));
p.addParameter('roiNames', struct(), @isstruct);

% Coordinate space option
p.addParameter('CoordSpace', 'RCS', @(x) any(strcmpi(x, {'LPH','RCS'})));

% ROI drawing style
p.addParameter('lw', 1, @(x) isnumeric(x) && isscalar(x));  % default line width
p.addParameter('fixedColor', [1 1 1], @(x) isnumeric(x) && numel(x)==3);  % default = white [R G B]

% UI Hiding controls (all default to TRUE)
p.addParameter('hideUI', true, @(x) islogical(x) || isnumeric(x));          
p.addParameter('hideCrosshairs', true, @(x) islogical(x) || isnumeric(x));  
p.addParameter('hideColorbar', true, @(x) islogical(x) || isnumeric(x));    
p.addParameter('hideROIText', true, @(x) islogical(x) || isnumeric(x));     
p.addParameter('hideButtons', true, @(x) islogical(x) || isnumeric(x));     
p.addParameter('hideVolName', true, @(x) islogical(x) || isnumeric(x));     

% For file path naming
p.addParameter('DesignmatName', '', @(x) ischar(x) || isstring(x));
p.addParameter('Modality', '',     @(x) ischar(x) || isstring(x));

p.parse(varargin{:});
opts = p.Results;

%% --------------------
% Setup
%% --------------------
scaleCoord = @(coord, scale) round(coord * scale);

if isempty(coords)
    coords = [99 99 129];
end
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Normalize fstem into a cell array
if ischar(fstem) || isstring(fstem)
    fstem = {char(fstem)};
end

% Get the showVol window
showVolWins = findall(groot, 'Type', 'figure', 'Tag', 'showVol');
if isempty(showVolWins)
    error('No showVol window found.');
end
myWin = showVolWins(1);
handles = guidata(myWin);

nVols = numel(handles.vols);

% Expand fstem if needed
if numel(fstem) == 1
    fstem = repmat(fstem, 1, nVols);
elseif numel(fstem) ~= nVols
    error('fstem must be a single string or match number of vols (%d).', nVols);
end

% Orientation state
orientation = handles.ORIENTATION;
set(0, 'CurrentFigure', myWin);
axes(handles.axes1);
hold(handles.axes1, 'on');

%% --------------------
% Hide UI elements (one-time global hiding)
%% --------------------
if opts.hideUI
    if opts.hideCrosshairs
        set(handles.toggleLines1, 'ForegroundColor', 'k');
        set(handles.axes1_X, 'Visible', 'off');
    end
    if opts.hideColorbar && isfield(handles, 'colorbar') && isgraphics(handles.colorbar)
        set(handles.colorbar, 'Visible', 'off');
        set(allchild(handles.colorbar), 'Visible', 'off');
    end
    if opts.hideButtons
        delete(findall(myWin,'Style','pushbutton','String','O'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines1'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines2'))
        delete(findall(myWin,'Style','togglebutton','Tag','toggleLines3'))
        if isfield(handles,'overlay_opacity_slider')
            delete(handles.overlay_opacity_slider);
        end
    end
end

updateDisplay_newvol(handles, true);

%% --------------------
% Loop through volumes & coords
%% --------------------
for v = 1:nVols
    fprintf('Processing vol %d/%d…\n', v, nVols);

    % Set current volume state
    set(0, 'CurrentFigure', myWin);
    tmp = find(handles.hideVol == 0);
    handles.currentVol = tmp(v);
    handles.currentVolPB = handles.volPB(handles.currentVol);
    set(handles.currentVolPB, 'ForegroundColor', 'r');
    handles.Mvxl2lph = handles.vols{handles.currentVol}.Mvxl2lph;

    handles.ORIENTATION = orientation;
    handles = displayNewOrientation(handles);  
    guidata(myWin, handles);                    
    updateDisplay_newvol(handles, true);   

    V = handles.vols{handles.currentVol};

    %% --------------------
    % Coord loop
    %% --------------------
    for ci = 1:size(coords,1)
        coord = coords(ci,:);    

        ud = struct();
        if strcmpi(opts.CoordSpace, 'LPH')
            ud.lph = coord;
        else
            Mvxl2lph = handles.Mvxl2lph_atlas;
            lph_coord = Mvxl2lph * [coord 1]';  % convert to LPH for showVol
            ud.lph = lph_coord(1:3)';
        end
        ud.M = handles.Mvxl2lph_atlas;
        ud.sender = myWin;
        set(myWin, 'UserData', ud);

        % Handle coordinate space conversion
        if strcmpi(opts.CoordSpace, 'LPH')
            % Convert LPH → voxel indices
            Mlph2vxl = inv(handles.Mvxl2lph_atlas);
            targetRCS = Mlph2vxl * [coord 1]';   
            handles.rr = targetRCS(1);
            handles.cc = targetRCS(2);
            handles.ss = targetRCS(3);
        else
            % Already in voxel indices
            handles.rr = coord(1);
            handles.cc = coord(2);
            handles.ss = coord(3);
        end

        % Refresh display for this coordinate
        handles = updateDisplay_newvol(handles, true);

        % Hide UI elements that need hiding every loop
        if opts.hideUI
            if opts.hideVolName && isfield(handles, 'volName_text') && isgraphics(handles.volName_text)
                set(handles.volName_text, 'Visible', 'off');
            end
            if opts.hideROIText
                roiFields = {'aseg','fiber','aparcaseg','thalamus','pauli'};
                for rf = 1:numel(roiFields)
                    if isfield(handles.anat, roiFields{rf}) && ...
                       isfield(handles.anat.(roiFields{rf}), 'roi_text')
                        try
                            set(handles.anat.(roiFields{rf}).roi_text, 'Visible', 'off');
                        catch
                        end
                    end
                end
            end
        end

        % ROI overlay logic
        if opts.drawROI && ~isempty(opts.parcellation)
            for a = 1:numel(opts.parcellation)
                atlas = opts.parcellation{a};
                if ~isfield(handles.anat, atlas)
                    warning('Parcellation "%s" not found in handles.anat. Skipping.', atlas);
                    continue;
                end

                atlasStruct = handles.anat.(atlas);

                if isfield(opts.roiNames, atlas)
                    roi_list = opts.roiNames.(atlas);
                else
                    roi_list = {};
                end
                if ischar(roi_list) || isstring(roi_list)
                    roi_list = {roi_list};
                end

                expanded_roi_list = {};
                for r = 1:numel(roi_list)
                    keyword = lower(roi_list{r});
                    switch keyword
                        case {'subcortical_all','all_fiber','all_aparcaseg','all_thalamus','all_pauli'}
                            expanded_roi_list = [expanded_roi_list; atlasStruct.uiNames(2:end)];
                        otherwise
                            expanded_roi_list = [expanded_roi_list; roi_list{r}];
                    end
                end

                keep_idx = find(ismember(atlasStruct.uiNames, expanded_roi_list));

                atlasStruct.uiRoiOverlaySelected = unique([atlasStruct.uiRoiOverlaySelected, keep_idx]);
                atlasStruct.uiRoiOverlayIdx = atlasStruct.uiRoiIdx(atlasStruct.uiRoiOverlaySelected);
                atlasStruct.uiRoiOverlayImg = atlasStruct.prob(:,:,:,atlasStruct.uiRoiOverlayIdx);

                if isfield(atlasStruct, 'params')
                    atlasStruct.params.showOverlay = true;
                    atlasStruct.params.showOutline = true;
                end

                handles.anat.(atlas) = atlasStruct;
            end

            handles.anat.params.showOverlay = true;
            handles.anat.params.showOutline = true;
            handles.anat.params.overlayAlpha = 0;
            set(handles.overlay_cb, 'Value', true);

            % Draw ROI outlines
            sc = calculateScale(handles);
            switch handles.ORIENTATION
                case 1 
                    sliceIndex = scaleCoord(handles.ss, sc.s);
                    ima = squeeze(V.imgs(:,:,sliceIndex));
                    rr = 1:size(ima,1); cc = 1:size(ima,2); ss = sliceIndex;
                    sliceLabel = 'CR';
                case 2
                    sliceIndex = scaleCoord(handles.cc, sc.c);
                    ima = squeeze(V.imgs(:,sliceIndex,:));
                    rr = 1:size(ima,1); ss = 1:size(ima,2); cc = sliceIndex;
                    sliceLabel = 'SR';
                case 3
                    sliceIndex = scaleCoord(handles.rr, sc.r);
                    ima = squeeze(V.imgs(sliceIndex,:,:));
                    cc = 1:size(ima,1); ss = 1:size(ima,2); rr = sliceIndex;
                    sliceLabel = 'CS';
            end
            [~, outline] = anatomyAddRoiOverlay_save_images(handles, ima, rr, cc, ss, 'outline');
            handles = anatomyDrawRoiOutline_save_images(handles, outline, sliceLabel, opts.lw, opts.fixedColor);
            guidata(myWin, handles);
        end

        drawnow; pause(0.3);
        
        % NOTE: Colorbar setup moved to RIGHT BEFORE screenshot to prevent override

        % % -------------------------
        % % Save screenshots
        % % -------------------------
        % if strcmpi(opts.CoordSpace, 'RCS')
        %     % 3-digit padded for RCS (e.g., 096)
        %     coord_str = sprintf('X%03d_Y%03d_Z%03d', coord(1), coord(2), coord(3));
        % else
        %     % 2-digit padded for LPH, with N for negatives (e.g., N03)
        %     coord_parts = cell(1,3);
        %     for ii = 1:3
        %         if coord(ii) < 0
        %             coord_parts{ii} = sprintf('N%02d', abs(coord(ii)));
        %         else
        %             coord_parts{ii} = sprintf('%02d', coord(ii));
        %         end
        %     end
        %     coord_str = sprintf('X%s_Y%s_Z%s', coord_parts{:});
        % end

        % switch handles.ORIENTATION
        %     case 1, orientation_name = 'coronal';
        %     case 2, orientation_name = 'sagittal';
        %     case 3, orientation_name = 'axial';
        %     otherwise, orientation_name = sprintf('orient%d', handles.ORIENTATION);
        % end

        % outdir_full = fullfile(outdir, orientation_name);
        % if ~exist(outdir_full, 'dir'), mkdir(outdir_full); end

        % vol_str = sprintf('vol%03d', v);

        % if ~isempty(fstem{v})
        %     % coord_fstem_vol###
        %     fname_base = sprintf('%s_%s_%s', coord_str, fstem{v}, vol_str);
        % else
        %     % coord_vol###
        %     fname_base = sprintf('%s_%s', coord_str, vol_str);
        % end

        % % Final PNG path
        % fname_out = fullfile(outdir_full, [fname_base '.png']);

        % shot = getframe(handles.axes1);
        % imwrite(shot.cdata, fname_out);
        % fprintf(' → Saved: %s\n', fname_out);

        % -------------------------
        % Save screenshots
        % -------------------------
        if strcmpi(opts.CoordSpace, 'RCS')
            coord_str = sprintf('X%03d_Y%03d_Z%03d', coord(1), coord(2), coord(3));
        else
            coord_parts = cell(1,3);
            for ii = 1:3
                if coord(ii) < 0, coord_parts{ii} = sprintf('N%02d', abs(coord(ii)));
                else,             coord_parts{ii} = sprintf('%02d',  coord(ii));
                end
            end
            coord_str = sprintf('X%s_Y%s_Z%s', coord_parts{:});
        end

        switch handles.ORIENTATION
            case 1, orientation_name = 'coronal';
            case 2, orientation_name = 'sagittal';
            case 3, orientation_name = 'axial';
            otherwise, orientation_name = sprintf('orient%d', handles.ORIENTATION);
        end


        % Pull names from opts (and give safe fallbacks)
        designmatName = char(opts.DesignmatName);
        if isempty(designmatName), designmatName = 'designmat_unknown'; end
        modalityName  = char(opts.Modality);
        if isempty(modalityName),  modalityName  = 'modality_unknown';  end

        % directory layout:
        % <screens_root>/screenshots/<designmat>/<modality>/<orientation>/<coords>/
        outdir_full = fullfile(outdir, 'screenshots_with_colorbar', designmatName, modalityName, orientation_name, coord_str);
        if ~exist(outdir_full, 'dir'), mkdir(outdir_full); end

        % Keep vol number
        vol_str = sprintf('vol%03d', v);

        % filename: <modality>_<coords>_<fstem>_<vol###>.png
        if ~isempty(fstem{v})
            fname_base = sprintf('%s_%s_%s_%s', modalityName, coord_str, fstem{v}, vol_str);
        else
            fname_base = sprintf('%s_%s_%s', modalityName, coord_str, vol_str);
        end

        fname_out = fullfile(outdir_full, [fname_base '.png']);

        % colorbar setup 
        if ~opts.hideColorbar
            fprintf('DEBUG: Starting colorbar creation for volume %d\n', handles.currentVol);
            
            % Get current volume
            V_current = handles.vols{handles.currentVol};
            
            % Debug volume properties
            if isfield(V_current, 'colormap')
                fprintf('DEBUG: Volume has colormap, size: %dx%d\n', size(V_current.colormap, 1), size(V_current.colormap, 2));
            else
                fprintf('DEBUG: Volume missing colormap field\n');
            end
            
            if isfield(V_current, 'limits')
                fprintf('DEBUG: Volume has limits: [%.2f, %.2f]\n', V_current.limits(1), V_current.limits(2));
            else
                fprintf('DEBUG: Volume missing limits field\n');
            end
            
            % Check if current volume has colormap
            if isfield(V_current, 'colormap') && isfield(V_current, 'limits')
                
                % Get main axes position for reference
                main_pos = get(handles.axes1, 'Position');
                fprintf('DEBUG: Main axes position: [%.3f %.3f %.3f %.3f]\n', main_pos);
                
                % Position colorbar on left side of brain, make it longer
                cb_pos = [0.05, 0.2, 0.025, 0.6];  % [x, y, width, height] - left side, longer
                fprintf('DEBUG: Creating colorbar at position: [%.3f %.3f %.3f %.3f]\n', cb_pos);
                
                % Delete any existing colorbar
                try
                    existing_cb = findobj(myWin, 'Tag', 'MyColorbar');
                    if ~isempty(existing_cb)
                        delete(existing_cb);
                        fprintf('DEBUG: Deleted existing colorbar\n');
                    end
                catch
                    fprintf('DEBUG: No existing colorbar to delete\n');
                end
                
                % Create new colorbar axes
                cb_axes = axes('Parent', myWin, ...
                              'Position', cb_pos, ...
                              'Tag', 'MyColorbar');
                fprintf('DEBUG: Created colorbar axes\n');
                
                % Get colormap
                cmap = V_current.colormap;
                if size(cmap, 2) > 3
                    cmap = cmap(:, 1:3);
                end
                nColors = size(cmap, 1);
                fprintf('DEBUG: Using colormap with %d colors\n', nColors);
                
                % Create simple colorbar data
                cb_data = repmat((1:nColors)', 1, 5);
                
                % Display the colorbar
                try
                    h_img = imagesc(cb_data, 'Parent', cb_axes);
                    fprintf('DEBUG: Created imagesc successfully\n');
                    
                    % Apply colormap to this specific axes
                    colormap(cb_axes, cmap);
                    fprintf('DEBUG: Applied colormap\n');
                    
                    % Configure axes
                    set(cb_axes, 'YDir', 'normal');
                    set(cb_axes, 'XTick', []);
                    set(cb_axes, 'YLim', [1 nColors]);
                    set(cb_axes, 'YAxisLocation', 'right');  % Numbers on RIGHT side of colorbar
                    
                    % Simple tick setup - ONLY min and max values
                    tick_pos = [1, nColors];  % Only bottom and top
                    tick_vals = [V_current.limits(1), V_current.limits(2)];  % Min and max values

                    % formatting: use up to 6 decimal places if needed
                    format_val = @(x) sprintf('%.6f', x);  % Start with 6 decimals
                    tick_labels = arrayfun(@(x) regexprep(format_val(x), '\.?0+$', ''), tick_vals, 'UniformOutput', false);
                    
                    set(cb_axes, 'YTick', tick_pos);
                    set(cb_axes, 'YTickLabel', tick_labels);
                    set(cb_axes, 'FontSize', 9);
                    set(cb_axes, 'FontWeight', 'bold');
                    set(cb_axes, 'Color', 'none');  % Transparent background
                    set(cb_axes, 'XColor', 'w');    % White axis color
                    set(cb_axes, 'YColor', 'w');    % White tick marks and labels
                    
                    fprintf('DEBUG: Set ticks: %s\n', mat2str(tick_pos));
                    fprintf('DEBUG: Set labels: %s\n', strjoin(tick_labels, ', '));
                    
                    % Make sure it's visible
                    set(cb_axes, 'Visible', 'on');
                    set(h_img, 'Visible', 'on');
                    
                    % Force to front
                    uistack(cb_axes, 'top');
                    
                    fprintf('DEBUG: Colorbar setup complete\n');
                    
                catch ME
                    fprintf('DEBUG: Error creating colorbar: %s\n', ME.message);
                end
                
                % Multiple rendering attempts
                fprintf('DEBUG: Starting rendering\n');
                drawnow;
                pause(0.3);
                drawnow;
                fprintf('DEBUG: Rendering complete\n');
                
            else
                fprintf('DEBUG: Volume missing required fields - cannot create colorbar\n');
            end
        else
            fprintf('DEBUG: Colorbar creation skipped (hideColorbar = true)\n');
        end

        % Capture just the main axes (brain image only)
        fprintf('DEBUG: Capturing screenshot of main axes only\n');
        shot = getframe(handles.axes1);
        imwrite(shot.cdata, fname_out);
        fprintf(' → Saved: %s\n', fname_out);
    end
end

end