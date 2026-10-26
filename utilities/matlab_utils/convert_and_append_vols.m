function vols = convert_and_append_vols(vols, vol_stat, measure, stat_name, colnames_model, index, limits, pval, cmap, bg_vol, CSFprob, method, atlasVersion, vol_display_mode, gray_limits)
%CONVERT_AND_APPEND_VOLS Upsample a 2mm stat map and append showVol entries.
%
%   vols = convert_and_append_vols(vols, vol_stat, measure, stat_name, ...
%       colnames_model, index, limits, pval, cmap, bg_vol, CSFprob, method, ...
%       atlasVersion, vol_display_mode, gray_limits)
%
%   vol_stat         - 2mm masked stat volume (3-D) in ABCD atlas geometry
%   measure          - imaging modality label (e.g. 'rni')
%   stat_name        - figure title fragment
%   colnames_model   - cell of predictor names (index selects label)
%   index            - predictor index for naming (scalar)
%   limits           - [min max] or [min max threshold] for color overlay
%   gray_limits      - optional [min max] for grayscale entry; default = limits(1:2)
%   vol_display_mode - 'both' | 'color' | 'gray' (default 'both')
%
%   See also: convertFEMAVols, process_wald_stats, process_spline_slices

if nargin < 1 || isempty(vols)
    vols = {};
end
if nargin < 14 || isempty(vol_display_mode)
    vol_display_mode = 'both';
end
if nargin < 13 || isempty(atlasVersion)
    atlasVersion = 'ABCD3';
end
if nargin < 12 || isempty(method)
    method = 'linear';
end
if nargin < 11
    CSFprob = [];
end
if nargin < 10 || isempty(bg_vol)
    bg_vol = atlas_T1(atlasVersion, true);
elseif ischar(bg_vol) || isstring(bg_vol)
    bg_vol = load(char(bg_vol));
end
if nargin < 9 || isempty(cmap)
    cmap = mmil_cmap_blueblackred;
end
if nargin < 8
    pval = false;
end
if nargin < 7 || isempty(limits)
    statmax = max(abs(vol_stat(isfinite(vol_stat(:)))));
    limits = [-statmax statmax];
end
if nargin < 15 || isempty(gray_limits)
    gray_limits = limits(1:2);
else
    gray_limits = gray_limits(1:2);
end
if nargin < 6 || isempty(index)
    index = 1;
end
if ischar(colnames_model) || isstring(colnames_model)
    colnames_model = {char(colnames_model)};
end

vol_display_mode = lower(char(vol_display_mode));
switch vol_display_mode
    case {'both', 'color', 'gray'}
    otherwise
        error('convert_and_append_vols: vol_display_mode must be both, color, or gray.');
end

template_vol = atlas_T1(atlasVersion);

switch method
    case 'linear'
        statimgs = upsample_volume(vol_stat);
    case 'nearest'
        statimgs = upsample_volume_nearest(vol_stat);
    otherwise
        error('convert_and_append_vols: unknown method ''%s''.', method);
end

if ~isempty(CSFprob)
    statimgs = maskCSF(statimgs, CSFprob, atlasVersion);
end
statimgs = maskBrain(statimgs, atlasVersion);

if numel(limits) >= 3
    thresh = limits(3);
    thresh_str = sprintf(', > %s', num2str(thresh, 2));
else
    thresh = [];
    thresh_str = '';
end

if iscell(stat_name)
    stat_name = stat_name{1};
end
if numel(colnames_model) >= index
    formula = sprintf('%s ~ %s', measure, colnames_model{index});
else
    formula = sprintf('%s ~ %s', measure, stat_name);
end

if pval
    stat_str = 'log p';
else
    stat_str = stat_name;
end

if any(strcmp(vol_display_mode, {'both', 'gray'}))
    vol_gray = template_vol;
    vol_gray.imgs = statimgs;
    vol_gray.name = sprintf('%s [%s]', formula, stat_name);
    vol_gray.statname = stat_name;
    vol_gray.limits = gray_limits;
    vol_gray.minI = gray_limits(1);
    vol_gray.maxI = gray_limits(2);
    vol_gray.colormap = gray(128);
    vols{end + 1} = vol_gray; %#ok<AGROW>
end

if any(strcmp(vol_display_mode, {'both', 'color'}))
    [vol_stat_rgb, cmap_out] = Get_volOverlay2(statimgs, limits(1:2), bg_vol, pval, cmap);
    vol_color = template_vol;
    vol_color.imgs = vol_stat_rgb.imgs;
    vol_color.limits = limits(1:2);
    if pval && any(isinf(vol_color.limits))
        vol_color.limits = [-322 322];
    end
    vol_color.colormap = cmap_out;
    vol_color.statname = stat_str;
    vol_color.name = sprintf('%s [%s +/- %s%s]', formula, stat_str, num2str(vol_color.limits(2), 3), thresh_str);
    vols{end + 1} = vol_color; %#ok<AGROW>
end

end
