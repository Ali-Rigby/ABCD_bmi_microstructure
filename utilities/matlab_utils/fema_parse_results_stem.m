function [datatype, fstem_imaging] = fema_parse_results_stem(filename)
%FEMA_PARSE_RESULTS_STEM Parse datatype and fstem_imaging from FEMA result filename.
%
% Supports:
%   FEMA_wrapper_output_voxel_rni.mat              (legacy, single file)
%   FEMA_wrapper_output_voxel_rni_estimates.mat    (FEMA 4.x + outPrefix)
%
% Returns empty strings if the name cannot be parsed (e.g. bare FEMA_estimates.mat).

[~, base, ~] = fileparts(filename);

if endsWith(base, '_estimates')
    base = extractBefore(base, '_estimates');
end

prefix = 'FEMA_wrapper_output_';
if ~startsWith(base, prefix)
    datatype = '';
    fstem_imaging = '';
    return;
end

middle = extractAfter(base, prefix);
parts = split(middle, '_');
if isempty(parts) || strlength(parts(1)) == 0
    datatype = '';
    fstem_imaging = '';
    return;
end

datatype = char(parts(1));
if numel(parts) > 1
    fstem_imaging = char(strjoin(parts(2:end), '_'));
else
    fstem_imaging = '';
end
end
