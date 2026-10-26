function minMax = pheno_basis_minmax(valvec, paths)
%PHENO_BASIS_MINMAX minMax argument for createBasisFunctions (s3), per pheno profile.
%
%   bmip         - fixed 0:100 BMI percentile grid (101 evaluation rows; step 1)
%   perc_median  - data-driven grid via pheno_eval_grid (BMI_IMAGING_PHENO_GRID_STEP)

switch paths.pheno_variable
    case 'bmip'
        minMax = [0 100 101];
    case 'perc_median'
        [~, minMax] = pheno_eval_grid(valvec);
    otherwise
        error('pheno_basis_minmax: unsupported BMI_PHENO_VARIABLE ''%s''.', paths.pheno_variable);
end

end
