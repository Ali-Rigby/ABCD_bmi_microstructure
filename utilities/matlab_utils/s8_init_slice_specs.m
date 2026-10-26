function slice_specs = s8_init_slice_specs(analysis_profile)
%S8_INIT_SLICE_SPECS Default slice_specs for s8 / showVol extended runs.
%
% Age spline slices:
%   Main analysis and ve11_nove11_* sensitivity profiles use 8.5, 9:17, 17.8
%   (same ages as original full-sample indices on 8.3:0.1:17.8).

slice_specs = struct();

slice_specs.bmi.spline = struct('mode', 'indices', 'indices', [2, 11:10:91, 100]);
% Full-sample / VE11 vs noVE11 matched arms
slice_specs.age.spline = struct('mode', 'values', 'values', [8.5, 9:17, 17.8]);
slice_specs.age.deriv  = struct('mode', 'values', 'values', [8.5, 9:17, 17.8]);

% Unused here but kept for callers that check the field
if nargin < 1
    analysis_profile = ''; %#ok<NASGU>
end

end
