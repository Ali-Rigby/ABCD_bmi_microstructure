function g = bmi_showvol_gray_limits(gray_spec, color_lims, xname)
%BMI_SHOWVOL_GRAY_LIMITS Resolve grayscale [min max] for showVol maps.
%
%   g = bmi_showvol_gray_limits([], color_lims)
%       Use color_lims for grayscale (default when no gray override is set).
%
%   g = bmi_showvol_gray_limits([0 8], color_lims)
%       Explicit grayscale override (e.g. phenolim_gray).
%
%   g = bmi_showvol_gray_limits(gwaldlims_dict, color_lims, xname)
%       Per-xname struct override; missing/empty field falls back to color_lims.

if nargin < 2 || isempty(color_lims)
    error('bmi_showvol_gray_limits: color_lims required.');
end
color_lims = color_lims(1:2);

if nargin < 1 || isempty(gray_spec)
    g = color_lims;
    return;
end

if isstruct(gray_spec)
    if nargin < 3 || isempty(xname)
        error('bmi_showvol_gray_limits: xname required when gray_spec is a struct.');
    end
    if isfield(gray_spec, xname) && ~isempty(gray_spec.(xname))
        g = gray_spec.(xname);
        g = g(1:2);
    else
        g = color_lims;
    end
else
    g = gray_spec(1:2);
end

end
