function [xvec, minMax] = pheno_eval_grid(valvec, step)
%PHENO_EVAL_GRID Clean integer x-grid for %from-median pheno maps/plots.
%
% Rounds data range outward to multiples of STEP so axis labels are whole numbers.
%   step = 5  ->  ..., -30, -25, -20, ..., 135, 140
%   step = 1  ->  every integer from floor(min) to ceil(max)
%
% Set BMI_IMAGING_PHENO_GRID_STEP before s3/s8/s9 for perc_median (ignored for bmip).
%
% Outputs:
%   xvec   - evaluation x values (row vector)
%   minMax - [minX maxX n] for createBasisFunctions minMax argument

if nargin < 2 || isempty(step)
    step = str2double(getenv('BMI_IMAGING_PHENO_GRID_STEP'));
    if isnan(step) || step <= 0
        step = str2double(getenv('BMI_MEDIAN_PHENO_GRID_STEP'));
    end
    if isnan(step) || step <= 0
        step = 5;
    end
end

valvec = valvec(:);
minX = step * floor(min(valvec) / step);
maxX = step * ceil(max(valvec) / step);
xvec = minX:step:maxX;
minMax = [minX, maxX, numel(xvec)];

end
