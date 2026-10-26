function process_voxel_data(fname_results, dirname_out, linear_terms, xnames, dirname_desmat, bf_dir, varargin)
%PROCESSVOXELDATA Wald/spline voxel stats using FEMA estimates + mask loader.
    p = inputParser;
    addParameter(p, 'Datatype', '', @ischar);
    addParameter(p, 'FstemImaging', '', @ischar);
    addParameter(p, 'DirnameImaging', '', @ischar);
    addParameter(p, 'BfDirShared', '', @ischar);
    parse(p, varargin{:});

    bf_dir_shared = p.Results.BfDirShared;
    if isempty(bf_dir_shared)
        bf_dir_shared = bf_dir;
    end

    R = fema_load_results('ResultsFile', fname_results, ...
        'Datatype', p.Results.Datatype, ...
        'FstemImaging', p.Results.FstemImaging, ...
        'DirnameImaging', p.Results.DirnameImaging);

    beta_hat = R.beta_hat;
    coeffCovar = R.coeffCovar;
    mask = R.mask;
    colnames_model = R.colnames_model;
    fname_results = R.fname_results;

    for i = 1:numel(linear_terms)
        linear_term = linear_terms{i};
        xname       = xnames{i};

        has_linear  = any(strcmp(colnames_model, linear_term));
        has_splines = any(contains(colnames_model, ['bf_', xname]));
        if ~(has_linear || has_splines), continue; end

        fprintf('Computing voxelwise stats for: %s\n', xname);
        compute_voxel_stats(beta_hat, coeffCovar, mask, colnames_model, dirname_out, ...
                          linear_term, xname, fname_results, dirname_desmat, bf_dir, bf_dir_shared);
    end
end
