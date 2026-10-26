function xvec = pheno_bf_xvec(tbl_bf_pheno, paths)
%PHENO_BF_XVEC Evaluation x-axis aligned with table_bf_* row count (s3/s7/s8/s9).
%
%   bmip         - 0:100 integer BMI percentiles (101 points)
%   perc_median  - pheno_eval_grid (BMI_IMAGING_PHENO_GRID_STEP from env)

n_bf = height(tbl_bf_pheno);

switch paths.pheno_variable
    case 'bmip'
        xvec = (0:100)';
    case 'perc_median'
        vec_file = fullfile(paths.datafiles, paths.pheno_vec_file);
        val_table = readtable(vec_file, 'Delimiter', '\t');
        valvec = val_table.(val_table.Properties.VariableNames{1});
        xvec = pheno_eval_grid(valvec);
    otherwise
        error('pheno_bf_xvec: unsupported BMI_PHENO_VARIABLE ''%s''.', paths.pheno_variable);
end

if n_bf ~= numel(xvec)
    error(['pheno_bf_xvec: BF table has %d rows but expected %d for %s. Re-run s3.'], ...
        n_bf, numel(xvec), paths.pheno_variable);
end

end
