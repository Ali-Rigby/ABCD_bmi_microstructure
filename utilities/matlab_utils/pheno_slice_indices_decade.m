function idx = pheno_slice_indices_decade(phenovec, decade_step)
%PHENO_SLICE_INDICES_DECADE Grid indices for spline showVol slices.
%
%   idx = pheno_slice_indices_decade(phenovec)
%
%   Returns indices into phenovec for:
%     1) the minimum grid value
%     2) every multiple of decade_step (default 10) from the next decade
%        at/above the minimum up through the last decade at/below the maximum
%     3) the maximum grid value
%
%   Example (min = -27, max = 134, step = 1 grid):
%     values -27, -20, -10, 0, 10, ..., 130, 134

if nargin < 2 || isempty(decade_step)
    decade_step = 10;
end

phenovec = phenovec(:);
n = numel(phenovec);
if n == 0
    idx = [];
    return;
end

vmin = phenovec(1);
vmax = phenovec(end);

first_decade = ceil(vmin / decade_step) * decade_step;
last_decade = floor(vmax / decade_step) * decade_step;
if first_decade <= last_decade
    decade_vals = first_decade:decade_step:last_decade;
else
    decade_vals = [];
end

target_vals = unique([vmin; decade_vals(:); vmax], 'stable');

idx = zeros(numel(target_vals), 1);
for i = 1:numel(target_vals)
    [~, idx(i)] = min(abs(phenovec - target_vals(i)));
end
idx = unique(idx);
idx = sort(idx);

% Guard against out-of-range (should not happen)
idx = idx(idx >= 1 & idx <= n);

end
